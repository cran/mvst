wQuantiles = function(x, weights, p){
 sort.list.x = sort.list(x)
 sort.x = x[sort.list.x]
 cum.weights = cumsum(weights[sort.list.x])
 q = numeric(length(p))
 for(ip in 1:length(p)) q[ip] = sort.x[which(cum.weights > p[ip])[1]]
 return(q)
}
