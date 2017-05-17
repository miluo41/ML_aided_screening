library(tidyverse);library(readxl);library(lubridate)

#read the kinetic data of cell fluorescence, induction OD, GC results
#header=T is default in read_xlsx function
#output in Hadley Wickham tidy format.

read_fluorescence<-function(file_name,signal_range,time_range){
        fluorsecence<-read_xlsx(file_name,range=signal_range)
        time<-read_xlsx(file_name,range=time_range,
                        col_types="date") %>%
                mutate(time=hour(.[[1]])*60+minute(.[[1]])) %>%
                select(time)
        output<-as_tibble(cbind(time['time'],fluorsecence)) %>%
                gather(key=well_num,value=reading,-time)
        output
}

read_induction<-function(file_name,induction_range){
        ind_OD<-as.matrix(read_xlsx(file_name,range=induction_range,
                                    col_types = 'numeric'))
        dim(ind_OD)<-c(96,1)
        plate<-outer(LETTERS[1:8],1:12,paste0)
        dim(plate)<-c(96,1)
        output<-as_tibble(cbind(plate,ind_OD)) %>%
                mutate(V2=as.numeric(V2)) %>%
                rename(well_num=V1,ind_OD=V2)
        output
}

read_GC<-function(file_name,well_range,GC_range){
        well<-read_xlsx(file_name,range=well_range)
        GC<-read_xlsx(file_name,range=GC_range,na='none')
        GC[is.na(GC)]<-0
        output<-as_tibble(cbind(well,GC))
        output

}