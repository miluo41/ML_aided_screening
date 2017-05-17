library(zoo);library(stringr)

extract_slope_516<-function(df,low=0.3,high=0.7){
        span<-max(df[['OD516']])-min(df[['OD516']])
        lower<-min(df[['OD516']])+0.2*span
        upper<-max(df[['OD516']])-0.2*span
        df_trim<-df[df[['OD516']] > lower & df[['OD516']] < upper,]
        coef(lm(OD516 ~ time,data=df_trim))[[2]]
}

extract_slope_545<-function(df,low=0.3,high=0.7){
        span<-max(df[['OD545']])-min(df[['OD545']])
        lower<-min(df[['OD545']])+0.2*span
        upper<-max(df[['OD545']])-0.2*span
        df_trim<-df[df[['OD545']] > lower & df[['OD545']] < upper,]
        coef(lm(OD545 ~ time,data=df_trim))[[2]]
}
take_auc<-function(df,y){
        a<-df[['time']]
        b<-df[[y]]
        id<-order(a)
        auc<-sum(diff(a[id])*rollmean(b[id],2))
        auc
}
