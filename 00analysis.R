#load files and packages
source('01import.R');source('02feature.R');
library(caret);library(e1071);library(corrplot);
library(leaps);library(ModelMetrics);library(tidyverse);
library(glmnet);library(ggthemes);
library(modelr)
#--------------------------------------------------------------------
#data wrangling
#fluorescence data
file_name<-"plate_1_fluorsecent.xlsx"
OD516_range<-"AZ51:CU191"
time_range<-"B51:B191"
OD516_df<-read_fluorescence(file_name,OD516_range,time_range)
names(OD516_df)[[3]]<-'OD516'
OD545_range<-"AZ200:CU340"
OD545_df<-read_fluorescence(file_name,OD545_range,time_range)
names(OD545_df)[[3]]<-'OD545'

fluo_df<-inner_join(OD516_df,OD545_df,by=c('time','well_num'))

# fluo_df %>% ggplot(aes(time,OD516))+
#         geom_line(aes(group=well_num,col=well_num))

#generate features
feature_df_1<-fluo_df %>%
        group_by(well_num) %>%
        nest() %>%
        mutate(slope.516=map(data,extract_slope_516),
               slope.545=map(data,extract_slope_545),
               int.516=map2(data,'OD516',take_auc),
               int.545=map2(data,'OD545',take_auc))%>%
        unnest(slope.516,slope.545,int.516,int.545,.drop=TRUE)

feature_df_2<-fluo_df %>%
        group_by(well_num) %>%
        summarize(max.516=max(OD516),
                  max.545=max(OD545),
                  min.516=min(OD516),
                  min.545=min(OD545),
                  median.516=median(OD516),
                  median.545=median(OD545))

#induction data
file_induction<-'plate_1_induction.xlsx'
induction_range<-'C26:N34'
induction_df<-read_induction(file_induction,induction_range) %>%
        filter(!str_sub(well_num,1,1) %in% c("A","B","C","D"))

#GC data
file_GC<-'plate_1_GC.xlsx'
well_range<-'G5:G53'
GC_range<-'I5:I53'
GC_df<-read_GC(file_GC,well_range,GC_range)
names(GC_df)<-c('well_num','peak.area')

#merge
main_df<-left_join(feature_df_1,feature_df_2,by='well_num') %>%
        left_join(.,induction_df,by='well_num') %>%
        left_join(.,GC_df,by='well_num')

#--------------------------------------------------------------------
# preprocess

# check for skewness
skew<-lapply(main_df[-1],skewness)

# #check for correlation
# 
corre.df<-main_df[,c(-1,-13)]
corre1<-cor(corre.df)
highCorr<-findCorrelation(corre1)
corre.df.low<-corre.df[,-highCorr]
corre2<-cor(corre.df.low)

corre.plot1<-corrplot::corrplot(corre1,order='hclust')
corre.plot2<-corrplot::corrplot(corre2,order='hclust')
df<-data.frame(corre.df.low,GC=main_df$peak.area)
#data splitting

# #data splitting
df.hist<-df %>% ggplot()+
        geom_histogram(aes(df$GC)) +
        xlab('GC Peak area (a.u.)')+
        geom_vline(xintercept=40)


set.seed(123)
inTrain <- createDataPartition(df$GC,p=0.8,list=FALSE)
training <- df[inTrain,]
testing <- df[-inTrain,]

#--------------------------------------------------------------------
#apply machine learning

#OLS
regfit.full<-regsubsets(GC~.,data=training)
reg.summary<-summary(regfit.full)
lm.fit<-train(GC~.,method='lm',data=training,
              preProc=c('BoxCox','scale','center'))
lm.pred<-predict(lm.fit,testing)
lm.rmse<-ModelMetrics::rmse(testing$GC,lm.pred)

#lasso
lasso.preprocess<-preProcess(training[,-6],
                             method=c('center','scale','BoxCox'))
training2<-predict(lasso.preprocess,training[,-6])
training2<-data.frame(training2,GC=training$GC)
testing2<-predict(lasso.preprocess,testing[,-6])
testing2<-data.frame(testing2,GC=testing$GC)
x<-model.matrix(GC~.,training2)[,-1]
y<-training2$GC
cv.out<-cv.glmnet(x,y,alpha=1)
lasso.fit<-glmnet(x,y,alpha=1,lambda = 10^seq(-2,10,length.out = 100))
test.x<-model.matrix(GC~.,testing2)[,-1]
lasso.pred<-predict(lasso.fit,s=cv.out$lambda.min,newx=test.x)
lasso.rmse<-ModelMetrics::rmse(testing$GC,lasso.pred)

#plsr
plsr.fit<-train(GC~.,training,method='pls',
                preProc=c('center','scale','BoxCox'),tuneLength=5)
plsr.pred<-predict(plsr.fit,testing)
plsr.rmse<-ModelMetrics::rmse(testing$GC,plsr.pred)
plsr.df<-data.frame(x=testing$GC,y=plsr.pred)
plsr.plot<-ggplot(plsr.df)+geom_point(aes(x,y),size=3)+theme_base()+
        xlab('observation (a.u.)')+ylab('rf.prediction (a.u.)')+
        ggtitle('rf Prediction')+ylim(0,80)+xlim(0,80)+
        coord_fixed()
#rf
set.seed(123)
ctrl=trainControl(method='repeatedcv',repeats = 2,number=10)
rf.fit<-train(GC~.,training,method='rf',trControl=ctrl,
              preProc=c('center','scale','BoxCox'),
              tuneLength=4)
rf.pred<-predict(rf.fit,testing)
rf.rmse<-ModelMetrics::rmse(testing$GC,rf.pred)
rf.df<-data.frame(x=testing$GC,y=rf.pred)
rf.plot<-ggplot(rf.df)+geom_point(aes(x,y),size=3)+theme_base()+
        xlab('observation (a.u.)')+ylab('rf.prediction (a.u.)')+
        ggtitle('rf Prediction')+ylim(0,80)+xlim(0,80)+
        coord_fixed()

#svm
set.seed(123)
svm.fit<-train(GC~.,training,
               method='svmRadial',
               preProc=c('center','scale'),
               tuneLength = 14,
               trControl = trainControl(method='cv'))
svm.pred<-predict(svm.fit,testing)
svm.rmse<-ModelMetrics::rmse(testing$GC,svm.pred)
svm.df=data.frame(x=testing$GC,y=svm.pred)
svm.plot<-ggplot(svm.df)+
        geom_abline(slope=1,intercept = 0,size=2,col='red3')+
        geom_point(aes(x,y),size=4)+theme_base()+
        xlab('observation (a.u.)')+ylab('svm.prediction (a.u.)')+
        ggtitle('SVM Prediction')+coord_fixed()+
        xlim(0,80)+ylim(0,80)+
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
svm.df2<-data.frame(x=rep(c(1:8),2),y=c(testing$GC,svm.pred),
                    z=rep(c('observation','prediction'),each=8))
svm.plot2<-ggplot(svm.df2)+
        geom_point(aes(x,y,group=z,col=z),size=4)+
        scale_color_manual(values=c('grey70','red4'))+
        theme_base()+
        xlab('sample ID')+ylab('Observation/Prediction (a.u.)')+
        ggtitle('SVM Prediction')+coord_fixed(ratio=0.1)+
        theme(legend.title=element_blank(),
              legend.text=element_text(size=15),
              legend.position=c(0.75,0.85),
              axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
#nnet
set.seed(123)
nnet.grid<-expand.grid(decay=c(0,0.01,0.1),
                       size=c(1:5),
                       bag=FALSE)
nnet.fit<-train(GC~.,training,
                method='avNNet',
                preProc=c('center','scale','pca'),
                tuneGrid=nnet.grid)
nnet.pred<-predict(nnet.fit,testing)
nnet.rmse<-ModelMetrics::rmse(testing$GC,nnet.pred)
#--------------------------------------------------------------------
# model comparison

model.df<-data.frame(model=c('lm','lasso','plsr','rf','svm','nnet'),
                     rmse=c(lm.rmse,lasso.rmse,
                     plsr.rmse,rf.rmse,svm.rmse,
                     nnet=nnet.rmse),
                     highlight=c(rep('non',4),'highlight','non'))
model.plot<-model.df %>% ggplot()+
        geom_col(aes(model,rmse,fill=highlight))+
        scale_fill_manual(values = c('red4','grey70'))+
        theme_base()+
        ggtitle('Model Comparison')+
        xlab('Model Name')+ylab('RMSE')+coord_fixed(ratio=0.15)+
        theme(legend.position = 'none')

#--------------------------------------------------------------------
#Bootstrapping SVM
complete.set<-rbind(training,testing)
index<-sample(48,80,replace = TRUE)
boot.set<-complete.set[index,] %>%
        add_predictions(svm.fit,'svm.pred') %>%
        add_residuals(svm.fit,'svm.resid')
boot.matrix<-matrix(boot.set$svm.resid,10,8)
boot.matrix2<-boot.matrix^2
boot.rmse<-apply(boot.matrix2,1,function(x) sqrt(sum(x)/length(x)))
svm.rmse.sd<-sd(boot.rmse)
svm.rmse.mean<-mean(boot.rmse)

#--------------------------------------------------------------------
png(file='corre1.png')
corrplot::corrplot(corre1,order='hclust')

png(file='corre2.png')
corrplot::corrplot(corre2,order='hclust')

png(file='svm_plot1.png')
svm.plot

png(file='svm_plot2.png')
svm.plot2

png(file='model_comparison.png')
model.plot
dev.off()
