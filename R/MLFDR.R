#' MLFDR: Mediation analysis using localFDR
#'
#' @param lfdr a vector of local false discovery rates, as output from the function localFDR
#' @param size Target False Discovery rate, default is 0.05
#'
#' @returns A vector of logicals indicating which hypotheses are rejected
#' @export
#'
#' @examples
MLFDR <- function(lfdr, size = 0.05){
  m = length(lfdr)
  st.lfdr<-sort(lfdr)
  k=1
  while(k<m && ((1/k)*sum(st.lfdr[1:k])) <= size){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<- lfdr<=lfdrk
  return(reject)
}
