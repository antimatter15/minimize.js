function minimize(X, f, length, varargin){
  // Minimize a differentiable multivariate function. 

  // Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )

  // where the starting point is given by "X" (D by 1), and the function named in
  // the string "f", must return a function value and a vector of partial
  // derivatives of f wrt X, the "length" gives the length of the run: if it is
  // positive, it gives the maximum number of line searches, if negative its
  // absolute gives the maximum allowed number of function evaluations. You can
  // (optionally) give "length" a second component, which will indicate the
  // reduction in function value to be expected in the first line-search (defaults
  // to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.

  // The function returns when either its length is up, or if no further progress
  // can be made (ie, we are at a (local) minimum, or so close that due to
  // numerical problems, we cannot get any closer). NOTE: If the function
  // terminates within a few iterations, it could be an indication that the
  // function values and derivatives are not consistent (ie, there may be a bug in
  // the implementation of your "f" function). The function returns the found
  // solution "X", a vector of function values "fX" indicating the progress made
  // and "i" the number of iterations (line searches or function evaluations,
  // depending on the sign of "length") used.

  // The Polack-Ribiere flavour of conjugate gradients is used to compute search
  // directions, and a line search using quadratic and cubic polynomial
  // approximations and the Wolfe-Powell stopping criteria is used together with
  // the slope ratio method for guessing initial step sizes. Additionally a bunch
  // of checks are made to make sure that exploration is taking place and that
  // extrapolation will not be unboundedly large.

  // See also: checkgrad 

  // Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).
  var SMALL = 10E-16
  var INT = 0.1;    // don't reevaluate within 0.1 of the limit of the current bracket
  var EXT = 3.0;                  // extrapolate maximum 3 times the current step-size
  var MAX = 20;                         // max 20 function evaluations per line search
  var RATIO = 10;                                       // maximum allowed slope ratio
  var SIG = 0.1; 
  var RHO = SIG/2; // SIG and RHO are the constants controlling the Wolfe-
  // Powell conditions. SIG is the maximum allowed absolute ratio between
  // previous and new slopes (derivatives in the search direction), thus setting
  // SIG to low (positive) values forces higher precision in the line-searches.
  // RHO is the minimum allowed fraction of the expected (from the slope at the
  // initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
  // Tuning of SIG (depending on the nature of the function to be optimized) may
  // speed up the minimization; it is probably not worth playing much with RHO.

  // The code falls naturally into 3 parts, after the initial line search is
  // started in the direction of steepest descent. 1) we first enter a while loop
  // which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
  // have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
  // enter the second loop which takes p2, p3 and p4 chooses the subinterval
  // containing a (local) minimum, and interpolates it, unil an acceptable point
  // is found (Wolfe-Powell conditions). Note, that points are always maintained
  // in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
  // conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
  // was a problem in the previous line-search. Return the best value so far, if
  // two consecutive line-searches fail, or whenever we run out of function
  // evaluations or line-searches. During extrapolation, the "f" function may fail
  // either with an error or returning Nan or Inf, and minimize should handle this
  // gracefully.

  // if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
  // if length>0, S='Linesearch'; else S='Function evaluation'; end 
  var red = 1;

  function dot(a, b){
    if(a.length != b.length || !(a.length > 0)){
      console.log(a, b)
      console.trace()
      throw "error dot";
    }
    var sum = 0;
    for(var i = 0; i < a.length; i++){
      sum += a[i] * b[i]
    }
    return sum
  }

  var i = 0;                                            // zero the run length counter
  var ls_failed = 0;                             // no previous line search has failed
  var feval = f(X, varargin);
  var f0 = feval[0], df0 = feval[1] // get function value and gradient

  var fX = [f0];
  i = i + (length<0);                                            // count epochs?!
  var s = df0.map(function(e){return -e}); 
  var d0 = -dot(s, s);           // initial search direction (steepest) and slope
  var x3 = red/(1-d0);                                  // initial step is red/(|s|+1)

  while(i < Math.abs(length)){                                      // while not finished
    i = i + (length>0);                                      // count iterations?!

    var X0 = X; 
    var F0 = f0; 
    var dF0 = df0;                   // make a copy of current values
    if(length>0){
      var M = MAX
    }else{
      var M = Math.min(MAX, -length-i)
    }

    while(true){                             // keep extrapolating as long as necessary
      var x2 = 0; 
      var f2 = f0; 
      var d2 = d0; 
      var f3 = f0; 
      var df3 = df0;
      var success = 0;
      // while ~success && M > 0
      while(!success && M > 0){
        try{
          M = M - 1; i = i + (length<0);                         // count epochs?!
          // [f3 df3] = feval(f, X+x3*s, varargin{:});
          // var feval = f(X + x3 * s, varargin); // TODO: implement scalar multiply
          var feval = f(X.map(function(e, i){ return e + x3 * s[i] }), varargin);
          f3 = feval[0]
          df3 = feval[1]
          if(isNaN(f3) || !isFinite(f3) || df3.some(function(e){return isNaN(e)}) || df3.some(function(e){return !isFinite(e)})){
            throw "error"
          }
          // if(isNaN(f3) || !isFinite(f3))
          // if isnan(f3) || isinf(f3) || any(isnan(df3)+isinf(df3)), error(''), end
          success = 1;
        }catch(err){                                // catch any error which occured in f
          x3 = (x2+x3)/2;                                  // bisect and try again
        }
      }
      if(f3 < F0){ // keep best values
        X0 = X.map(function(e, i){ return e + x3 * s[i] });
        F0 = f3; 
        dF0 = df3;
      }         
      var d3 = dot(df3, s)
      // d3 = df3'*s;                                                    // new slope
      if(d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0) break;  // are we done extrapolating?
        
      var x1 = x2; 
      var f1 = f2; 
      var d1 = d2;                        // move point 2 to point 1
      x2 = x3; f2 = f3; d2 = d3;                        // move point 3 to point 2
      var A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);                 // make cubic extrapolation
      var B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
      // x3 = x1-d1*(x2-x1)^2/(B+Math.sqrt(B*B-A*d1*(x2-x1))); // num. error possible, ok!
      // TODO: check to make sure that this is positive
      // console.log('---not---neg----', B * B - A * d1 * (x2 - x1))
      var Z = B + Math.sqrt(B * B - A * d1 * (x2 - x1))
      if(Z != 0){
        x3 = x1 - d1 * Math.pow(x2 - x1, 2) / Z
      }else{
        x3 = Infinity
      }

      // if ~isreal(x3) || isnan(x3) || isinf(x3) || x3 < 0 // num prob | wrong sign?
      if(!isFinite(x3) || isNaN(x3) || x3 < 0){
        x3 = x2*EXT;                                 // extrapolate maximum amount
      }else if(x3 > x2*EXT){                  // new point beyond extrapolation limit?
        x3 = x2*EXT;                                 // extrapolate maximum amount
      }else if(x3 < x2+INT*(x2-x1)){         // new point too close to previous point?
        x3 = x2+INT*(x2-x1);
      }
    }                                                       // end extrapolation
    while((Math.abs(d3) > -SIG * d0 || f3 > f0 * x3 * RHO * d0) && M > 0){

    // while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  // keep interpolating
      if(d3 > 0 || f3 > f0+x3*RHO*d0){                         // choose subinterval
        var x4 = x3; 
        var f4 = f3; 
        var d4 = d3;                      // move point 3 to point 4
      }else{
        var x2 = x3; 
        var f2 = f3; 
        var d2 = d3;                      // move point 3 to point 2
      }
      if(f4 > f0){
        x3 = x2-(0.5*d2*Math.pow(x4-x2, 2))/(f4-f2-d2*(x4-x2));  // quadratic interpolation
      } else {
        A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    // cubic interpolation
        B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
        x3 = x2+(Math.sqrt(B*B-A*d2*Math.pow(x4-x2, 2))-B)/A;        // num. error possible, ok!
      }
      if(isNaN(x3) || !isFinite(x3)){
        x3 = (x2+x4)/2;               // if we had a numerical problem then bisect
      }
      x3 = Math.max(Math.min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  // don't accept too close
      
      var feval = f(X.map(function(e, i){ return e + x3 * s[i] }), varargin);
      f3 = feval[0]
      df3 = feval[1]

      if(f3 < F0){ // keep best values
        X0 = X.map(function(e, i){ return e + x3 * s[i] });
        F0 = f3; 
        dF0 = df3;
      }         
      M = M - 1; i = i + (length<0);                             // count epochs?!
      d3 = dot(df3, s)
    }                                                     // end interpolation

    if(Math.abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0){          // if line search succeeded
      X = X.map(function(e, i){ return e + x3 * s[i] });


      f0 = f3; 
      fX.push(f0) //fX = [fX' f0]';                     // update variables
      // fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
      
      console.log(i, f0)

      var s_mul = (dot(df3, df3) - dot(df0, df3)) / dot(df0, df0)
      s = s.map(function(e, i){ return e * s_mul - df3[i] })
      df0 = df3;                                               // swap derivatives
      d3 = d0; 
      d0 = dot(df0, s);
      if(d0 > 0){                                      // new slope must be negative
        s = df0.map(function(e){return -e}); 

        d0 = -dot(s, s);                  // otherwise use steepest direction
      }
      x3 = x3 * Math.min(RATIO, d3/(d0-SMALL));          // slope ratio but max RATIO
      ls_failed = 0;                              // this line search did not fail
    }else{
      X = X0; f0 = F0; df0 = dF0;                     // restore best point so far
      if(ls_failed || i > Math.abs(length)){         // line search failed twice in a row
        break;                             // or we ran out of time, so we give up
      }
      s = df0.map(function(e){return -e});
      d0 = -dot(s, s);                                        // try steepest
      x3 = 1/(1-d0);                     
      ls_failed = 1;                                    // this line search failed
    }
  }
  return [X, fX, i]
}




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



function poisson_mix(theta, S){
  if(theta.length != 3) throw "yo theta gotta be 3";
  var ll1 = theta[0],
      ll2 = theta[1],
      lc = theta[2];
  S = S.map(function(e){return Math.max(1, Math.round(e))})
  // console.log(theta, S)
  var num = S.length
  var df = [0, 0, 0];
  if(!isFinite(ll1)) throw "walp lambda1 = 0";
  if(!isFinite(ll2)) throw "walp lambda2 = 0";
  var sumS = S.reduce(function(a,b){return a + b})
  // var c = lc.map(function(e){ return Math.exp(e) })
  var c = Math.exp(lc)
  var ecs = S.map(function(e){return Math.exp(c - e)})
  // var ecs = c.map(function(e, i){ return e - S[i] })
  var ex = S.map(function(e, i){
    return Math.exp(e - c - Math.exp(ll1) + e * ll1 + Math.exp(ll2) - e * ll2)
  })
  var ex1 = ex.map(function(e){return e + 1})
  function sum(a){
    var s = 0;
    for (var i = a.length - 1; i >= 0; i--) s += a[i];
    return s
  }

  // console.log('pecs', ex, S, ex1)
  var sumsigex = sum(ex.map(function(e, i){return e / ex1[i]}))
  var sumsigSex = sum(S.map(function(e, i){return e * ex[i] / ex1[i]}))
  var sumsigcex = sum(ex.map(function(e, i){return e * c / ex1[i]}))

  var logfactS = S.map(function(e){return 0})
  for(var i = 0; i < S.length; i++){
    for(var j = 0; j < S[i]; j++){
      logfactS[i] += Math.log(j + 1)
    }
  }

  sumlogfactS = sum(logfactS)
  // console.log(sumlogfactS, logfactS)
  // console.log('ecs', ecs, sumsigex, sumsigSex, sumsigcex)
  // console.log(sumlogfactS, sumsigex, sumsigSex, sumsigcex, 'sex')
  // console.log(sumlogfactS , ecs )
  
  var val = sumlogfactS + sum(ecs.map(function(e){return Math.log(e + 1)})) 
   - sum(S.map(function(e){
    return Math.log(Math.exp(e * ll1 - Math.exp(ll1)) + Math.exp(
        c - e + e * ll2 - Math.exp(ll2)
      ))
  }))

  var df1 = Math.exp(ll1) * sumsigex - sumsigSex
  var df2 = num * Math.exp(ll2) - sumS - Math.exp(ll2) * sumsigex + sumsigSex;
  
  var df3 = sum(ecs.map(function(e, i){
    return c * ecs[i] / (1 + ecs[i])
  })) - sum(S.map(function(e, i){
    return (c * Math.exp(c - e + e * ll2 - Math.exp(ll2))) / (
      Math.exp(e * ll1 - Math.exp(ll1)) + Math.exp(c - e + e * ll2 - Math.exp(ll2))
      )
  }))
  // console.log(val, df1, df2, df3)
  return [val, [df1, df2, df3]]
}



