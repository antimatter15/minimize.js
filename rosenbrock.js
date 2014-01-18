function rosenbrock(x){
  var a = x.slice(1)
  var b = x.slice(0, -1).map(function(e){ return Math.pow(e, 2) })
  var c = a.map(function(e, i){ return 100 * Math.pow(e - b[i], 2) })
  var d = x.slice(0, -1).map(function(e){ return Math.pow(1 - e, 2) })
  var f = 0;
  for(var i = 0; i < c.length; i++){
    f += c[i] + d[i]
  }
  var df = x.map(function(e){ return 0 })
  for(var i = 0; i < df.length - 1; i++){
    df[i] = -400 * x[i] * (x[i + 1] - Math.pow(x[i], 2)) - 2 * (1 - x[i])
  }
  for(var i = 1; i < df.length; i++){
    df[i] = df[i] + 200 * (x[i] - Math.pow(x[i - 1], 2))
  }
  return [f, df]
}