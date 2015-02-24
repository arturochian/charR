#' CharAnalysis
#'
#' Implements the CharAnalysis program to identify potential peaks in the charcoal record.
char <- function(z,smoothingtype,width,transform, peaksignal, threshold, globalusert, thresholdtype, local.window, percentile, pvalue){
  ## File input
  char.d <- data.frame(z$time, z$interpolated)
  colnames(char.d) <- c("years","series")

  ## Log transform if requested
  if (transform=="log"){
    char.d$interpolated <- log(char.d$series)
    char.d$interpolated[char.d$interpolated<0] <-0}
  if (transform=="identity"){
    char.d$interpolated <- char.d$series}

  ## smoothing algorithms to get the background

  if(smoothingtype=="lowess"){
    char.d$background<-lowess(x=1:length(char.d$interpolated),y=char.d$interpolated,f=width/length(char.d$interpolated),delta=width/2)$y}

  if(smoothingtype=="moving-average"){
    moving_average<-function(width){
      moving.average<-rep(1,width)/width
      return(moving.average)}
    moving.average<-moving_average(width)
    char.d$background<-filter(char.d$interpolated, moving.average, sides=2)}

  if(smoothingtype=="Savitsky-Golay"){
    savistsky_golay<-function(width){
      x<-1:width-width/2
      y<-max(x^2)-x^2
      sg<-y/sum(y)
      return(sg)}
    sg<-savistsky_golay(width)
    char.d$background<-filter(char.d$interpolated, sg, sides=2)}

  if(smoothingtype=="moving-median"){
    char.d$background<-med.filter(char.d$interpolated, width)$y}

  ## Peak signal:  either interpolated-background or interpolated/background
  if(peaksignal=="residual"){    char.d$peak <- char.d$interpolated-char.d$background }

  if(peaksignal=="index"){    char.d$peak <- char.d$interpolated/char.d$background }

  ### Flensed data

  char.d <- char.d[complete.cases(char.d),]
  attach(char.d)

  ### Threshold Shennanigans
  ## Global
  if(thresholdtype=="global"){
    if(threshold=="mixmodel"){
      mixmdl <- summary(Mclust(peak))
      t.temp <-qnorm(percentile, mean=mixmdl$mean[1], sd=sqrt(mixmdl$variance[1]))}

    if(threshold=="normal"){
      if(peaksignal=="residual"){
        mixmdl <- summary(Mclust(peak))
        t.temp <-qnorm(percentile, mean=0, sd=sqrt(mixmdl$variance[1]))}

      if(peaksignal=="index"){
        mixmdl <- summary(Mclust(peak))
        t.temp <-qnorm(percentile, mean=1, sd=sqrt(mixmdl$variance[1]))}}

    if(threshold=="specify"){
      t.temp <- globalusert}

    ## Thresholds after global
    t <- rep(NA, length(peak))
    t <- rep(t.temp, length(peak))}

  ## Local
  if(thresholdtype=="local"){
    t.temp <- c(rep(NA, length(peak)))
    if(threshold=="mixmodel"){

      for (i in 1:1){

        local.record <- peak[i:local.window]
        mixmdl <- Mclust(local.record)
        mu=summary(mixmdl)$mean[1]
        sigma=sqrt(summary(mixmdl)$variance[1])
        t.temp[i] <-qnorm(percentile, mu, sigma)}

      for (i in 2:ceiling(local.window/2)){

        local.record <- peak[i:local.window]
        mixmdl <- Mclust(local.record)
        mu=summary(mixmdl)$mean[1]
        sigma=sqrt(summary(mixmdl)$variance[1])
        t.temp[i] <-qnorm(percentile, mu, sigma)}

      for (i in (ceiling(local.window/2)+1):(length(t.temp)-ceiling(local.window/2))){

        local.record <- peak[ceiling(i-local.window/2):ceiling(i+local.window/2)]
        mixmdl <- Mclust(local.record)
        mu=summary(mixmdl)$mean[1]
        sigma=sqrt(summary(mixmdl)$variance[1])
        t.temp[i] <-qnorm(percentile, mu, sigma)}

      for (i in (length(t.temp)-ceiling(local.window/2)):length(t.temp)){
        local.record <- peak[(i-local.window):i]
        mixmdl <- Mclust(local.record)
        mu=summary(mixmdl)$mean[1]
        sigma=sqrt(summary(mixmdl)$variance[1])
        t.temp[i] <-qnorm(percentile, mu, sigma)}}

    if(threshold=="normal"){
      if(peaksignal=="residual"){

        for (i in 1:1){
          local.record <- peak[i:local.window]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 0, sigma)}

        for (i in 2:ceiling(local.window/2)){

          local.record <- peak[i:local.window]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 0, sigma)}

        for (i in (ceiling(local.window/2)+1):(length(t.temp)-ceiling(local.window/2))){

          local.record <- peak[ceiling(i-local.window/2):ceiling(i+local.window/2)]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 0, sigma)}

        for (i in (length(t.temp)-ceiling(local.window/2)):length(t.temp)){
          local.record <- peak[(i-local.window):i]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile,0, sigma)}}

      if(peaksignal=="index"){

        for (i in 1:1){
          local.record <- peak[i:local.window]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 1, sigma)}

        for (i in 2:ceiling(local.window/2)){

          local.record <- peak[i:local.window]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 1, sigma)}

        for (i in (ceiling(local.window/2)+1):(length(t.temp)-ceiling(local.window/2))){

          local.record <- peak[ceiling(i-local.window/2):ceiling(i+local.window/2)]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile, 1, sigma)}

        for (i in (length(t.temp)-ceiling(local.window/2)):length(t.temp)){
          local.record <- peak[(i-local.window):i]
          mixmdl <- Mclust(local.record)
          sigma=sqrt(summary(mixmdl)$variance[1])
          t.temp[i] <-qnorm(percentile,1, sigma)}
      }}

    ### Smoothing algorithm for local threshold same width as background smoother
    t <- rep(NA, length(peak))

    if(smoothingtype=="lowess"){
      t<-lowess(x=1:length(t.temp),y=t.temp,f=width/length(t.temp),delta=width/2)$y}

    if(smoothingtype=="moving-average"){
      moving_average<-function(width){
        moving.average<-rep(1,width)/width
        return(moving.average)}
      moving.average<-moving_average(width)
      t<-filter(t.temp, moving.average, sides=2)}

    if(smoothingtype=="Savitsky-Golay"){
      # Savitsky-Golay smoothing function
      savistsky_golay<-function(width){
        x<-1:width-width/2
        y<-max(x^2)-x^2
        sg<-y/sum(y)
        return(sg)}
      sg<-savistsky_golay(width)
      t<-filter(t.temp, sg, sides=2)}

    if(smoothingtype=="moving-median"){
      t<-med.filter(t.temp, width)$y}

  }

  ### Potential Peaks-everything that exceeded whatever threshold

  #   fireids <- c(rep(NA,length(years)))
  #
  #   for (i in 1:length(years)){
  #     fireids[i] <- ifelse(peak[i]>t[i], years[i], NA)}
  fireids <- ifelse(peak > t, years, NA)

  ### Screen and remove insignificant peaks via minimum count statistic

  d <- c(rep(NA, length(peak)-1))
  for (i in 1:length(d)){
    d[i]<-(interpolated[i]-(interpolated[i]+interpolated[i+1])*
             (0.5) - 0.5)/(sqrt((interpolated[i]+interpolated[i+1])*.25))
  }

  sig.fireids <- c(rep(NA,length(years)))

  for (i in 1:length(years)){
    sig.fireids[i] <- ifelse(d[i]>qnorm(pvalue), years[i], NA)}

  ### Shit to get bumped out
  return(list(
    years=years, background=background,
    interpolated=interpolated, peak=peak,
    fireids=fireids,
    sig.fireids=sig.fireids, t=t))
}
