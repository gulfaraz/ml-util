var mathjs = require("mathjs");

mathjs.config({
    number: "BigNumber",    // Default type of number:
                            // 'number' (default), 'BigNumber', or 'Fraction'
    precision: 64           // Number of significant digits for BigNumbers
});

mathjs.isComplex = function (value) {
    return value != null && value.isComplex;
};

mathjs.complexLarger = function (a, b) {
    var isLarger = true;
    if(mathjs.isComplex(a) || mathjs.isComplex(b)) {
        isLarger = mathjs.abs(a) > mathjs.abs(b) || (mathjs.abs(a) == mathjs.abs(b) && mathjs.arg(b) > mathjs.arg(b));
    } else {
        isLarger = mathjs.larger(a, b);
    }
    return isLarger;
};


mathjs.isFinite = function (value) {

  // check if value is Infinity
  if (value == Number.POSITIVE_INFINITY || value == Number.NEGATIVE_INFINITY){
    return false;
  }

  // check if any slot in the array is Infinity
  if (value instanceof Array){
    for (var i=0; i<value.length; i++){
      var item=value[i];
      if (item == Number.POSITIVE_INFINITY || item == Number.NEGATIVE_INFINITY || math.equals(item, Number.POSITIVE_INFINITY) || math.equals(item, Number.NEGATIVE_INFINITY)){
        return false;
      }
      return true;
    }
  }

  // unknown
  return true;
};

module.exports = function (f, X, options, P1, P2, P3, P4, P5) {
    console.log("Starting FMIN Conjugate Gradient Function...");
/*
function [X, fX, i] = fmincg(f, X, options, P1, P2, P3, P4, P5)
% Minimize a continuous differentialble multivariate function. Starting point
% is given by "X" (D by 1), and the function named in the string "f", must
% return a function value and a vector of partial derivatives. The Polack-
% Ribiere flavour of conjugate gradients is used to compute search directions,
% and a line search using quadratic and cubic polynomial approximations and the
% Wolfe-Powell stopping criteria is used together with the slope ratio method
% for guessing initial step sizes. Additionally a bunch of checks are made to
% make sure that exploration is taking place and that extrapolation will not
% be unboundedly large. The "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The function returns when either its length is up, or if no further
% progress can be made (ie, we are at a minimum, or so close that due to
% numerical problems, we cannot get any closer). If the function terminates
% within a few iterations, it could be an indication that the function value
% and derivatives are not consistent (ie, there may be a bug in the
% implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% Usage: [X, fX, i] = fmincg(f, X, options, P1, P2, P3, P4, P5)
%
% See also: checkgrad 
%
% Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
%
%
% (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
% 
% Permission is granted for anyone to copy, use, or modify these
% programs and accompanying documents for purposes of research or
% education, provided this copyright notice is retained, and note is
% made of any changes that have been made.
% 
% These programs and documents are distributed without any warranty,
% express or implied.  As the programs were written for research
% purposes only, they have not been tested to the degree that would be
% advisable in any important application.  All use of these programs is
% entirely at the user's own risk.
%
% [ml-class] Changes Made:
% 1) Function name and argument specifications
% 2) Output display
%
*/

    var length = 100; //else    length = 100;

    if(options && !isNaN(options.maxIterations)) { //if exist('options', 'var') && ~isempty(options) && isfield(options, 'MaxIter')
        length = options.maxIterations; //    length = options.MaxIter;
    } //end

    var CONSTANT = {
        "RHO" : 0.01, //RHO = 0.01;                            % a bunch of constants for line searches
        "SIG" : 0.5, //SIG = 0.5;       % RHO and SIG are the constants in the Wolfe-Powell conditions
        "INT" : 0.1, //INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
        "EXT" : 3.0, //EXT = 3.0;                    % extrapolate maximum 3 times the current bracket
        "MAX" : 20, //MAX = 20;                         % max 20 function evaluations per line search
        "RATIO" : 100 //RATIO = 100;                                      % maximum allowed slope ratio
    };

    var argumentList = Array.prototype.slice.call(arguments); //  argstr = [argstr, ',P', int2str(i)]; end argstr = [argstr, ')'];
    argumentList.splice(0, 3); //for i = 1:(nargin - 3)

    var updatedF = function(X) {
        var updatedArgumentList = argumentList.slice(0);
        updatedArgumentList.unshift(X);
        return f.apply(null, updatedArgumentList); //argstr = ['feval(f, X'];                      % compose string used to call function
    };

    var red = 1; // else red=1;
    var S = "Iteration - "; //S=['Iteration '];
    if(length && length.length === 2) {
        length = length[1];
        red = length[2]; //if max(size(length)) == 2, red=length(2); length=length(1);
    } // end
    var i = 0; //i = 0;                                            % zero the run length counter
    var lsFailed = 0; //ls_failed = 0;                             % no previous line search has failed
    var fX = []; //fX = [];
    var [ f1, df1 ] = updatedF(X); //[f1 df1] = eval(argstr);                      % get function value and gradient
    i = i + ((length<0) ? 1 : 0); //i = i + (length<0);                                            % count epochs?!
    var s = mathjs.multiply(-1, df1); //s = -df1;                                        % search direction is steepest
    var d1 = mathjs.multiply(mathjs.multiply(-1, mathjs.transpose(s)), s); //d1 = -s'*s;                                                 % this is the slope
    var z1 = mathjs.dotDivide(red, mathjs.subtract(1, d1)); //z1 = red/(1-d1);                                  % initial step is red/(|s|+1)
    while(i < mathjs.abs(length)) {//while i < abs(length)                                      % while not finished
        i = i + ((length>0) ? 1 : 0); //i = i + (length>0);                                      % count iterations?!
        var X0 = mathjs.clone(X); //X0 = X;
        var f0 = mathjs.clone(f1); // f0 = f1;
        var df0 = mathjs.clone(df1); // df0 = df1;                   % make a copy of current values
        var z1ProductS = mathjs.multiply(z1, s);
        z1ProductS = mathjs.reshape(z1ProductS, [ z1ProductS.length, 1 ]);
        X = mathjs.add(X, z1ProductS); //X = X + z1*s;                                             % begin line search
        var [ f2, df2 ] = updatedF(X); //[f2 df2] = eval(argstr);
        i = i + ((length<0) ? 1 : 0); //i = i + (length<0);                                          % count epochs?!
        var d2 = mathjs.multiply(mathjs.transpose(df2), s); //d2 = df2'*s;
        var f3 = mathjs.clone(f1); //f3 = f1;
        var d3 = mathjs.clone(d1); // d3 = d1;
        var z3 = mathjs.multiply(-1, z1); // z3 = -z1;             % initialize point 3 equal to point 1
        var M = mathjs.min(CONSTANT.MAX, mathjs.subtract(mathjs.multiply(-1, length), i)); // else M = min(MAX, -length-i);
        if(length > 0) { //if length>0
            M = CONSTANT.MAX; //, M = MAX;
        } // end
        var success = 0;
        var limit = -1; //  success = 0; limit = -1;                     % initialize quanteties
        while(1) { //while 1
            while( ( (mathjs.larger(f2, mathjs.add(f1, mathjs.multiply(z1, CONSTANT.RHO, d1)))) || (mathjs.larger(d2, mathjs.multiply(-1, CONSTANT.SIG, d1))) ) && (M > 0) ) {//while ((f2 > f1+z1*RHO*d1) || (d2 > -SIG*d1)) && (M > 0)
                limit = mathjs.clone(z1); //limit = z1;                                         % tighten the bracket
                var z2;
                if(mathjs.larger(f2, f1)) { //if f2 > f1
                    z2 = mathjs.subtract(z3, mathjs.divide(mathjs.multiply(0.5, d3, z3, z3), mathjs.add(mathjs.subtract(mathjs.multiply(d3, z3), f3), f2))); //z2 = z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3);                 % quadratic fit
                } else { //else
                    var A = mathjs.add(mathjs.multiply(6,mathjs.divide(mathjs.subtract(f2,f3),z3)),mathjs.multiply(3,mathjs.add(d2,d3))); //A = ((6*((f2-f3)/z3))+(3*(d2+d3)));                                 % cubic fit
                    var B = mathjs.subtract(mathjs.multiply(3, mathjs.subtract(f3, f2)), mathjs.multiply(z3, mathjs.add(d3, mathjs.multiply(2, d2)))); //B = 3*(f3-f2)-z3*(d3+2*d2);
                    z2 = mathjs.divide(mathjs.subtract(mathjs.sqrt(mathjs.subtract(mathjs.multiply(B, B), mathjs.multiply(A, d2, z3, z3))), B), A); //z2 = (sqrt(B*B-A*d2*z3*z3)-B)/A;       % numerical error possible - ok!
                } //end
                if(mathjs.isNaN(z2) || !mathjs.isFinite(z2)) { //if isnan(z2) || isinf(z2)
                    z2 = mathjs.divide(z3, 2); //z2 = z3/2;                  % if we had a numerical problem then bisect
                } //end
                z2 = mathjs.max(mathjs.min(z2, mathjs.multiply(CONSTANT.INT, z3)), mathjs.multiply(mathjs.subtract(1, CONSTANT.INT), z3)); //z2 = max(min(z2, INT*z3),(1-INT)*z3);  % don't accept too close to limits
                z1 = mathjs.add(z1, z2); //z1 = z1 + z2;                                           % update the step
                var z2ProductS = mathjs.multiply(z2, s);
                z2ProductS = mathjs.reshape(z2ProductS, [ z2ProductS.length, 1 ]);
                X = mathjs.add(X, z2ProductS); //X = X + z2*s;
                [ f2, df2 ] = updatedF(X); //[f2 df2] = eval(argstr);
                M = M - 1; //M = M - 1;
                i = i + ((length<0) ? 1 : 0); //i = i + (length<0);                           % count epochs?!
                d2 = mathjs.multiply(mathjs.transpose(df2), s); //d2 = df2'*s;
                z3 = mathjs.subtract(z3, z2); //z3 = z3-z2;                    % z3 is now relative to the location of z2
            } //end
            if(mathjs.larger(f2, mathjs.add(f1, mathjs.multiply(z1, CONSTANT.RHO, d1))) || mathjs.larger(d2, mathjs.multiply(-1, CONSTANT.SIG, d1))) { //if f2 > f1+z1*RHO*d1 || d2 > -SIG*d1
                break; //break;                                                % this is a failure
            } else if(mathjs.larger(d2, mathjs.multiply(CONSTANT.SIG, d1))) { //elseif d2 > SIG*d1
                success = 1;
                break; //success = 1; break;                                             % success
            } else if(M === 0) { //elseif M == 0
                break; //break;                                                          % failure
            } //end
            var z2;
            var A = mathjs.add(mathjs.multiply(6,mathjs.divide(mathjs.subtract(f2,f3),z3)),mathjs.multiply(3,mathjs.add(d2,d3))); //A = ((6*((f2-f3)/z3))+(3*(d2+d3)));                      % make cubic extrapolation
            var B = mathjs.subtract(mathjs.multiply(3, mathjs.subtract(f3, f2)), mathjs.multiply(z3, mathjs.add(d3, mathjs.multiply(2, d2)))); //B = 3*(f3-f2)-z3*(d3+2*d2);
            z2 = mathjs.divide(mathjs.multiply(-1, d2, z3, z3), mathjs.add(B, mathjs.sqrt(mathjs.subtract(mathjs.multiply(B, B), mathjs.multiply(A, d2, z3, z3))))); //z2 = -d2*z3*z3/(B+sqrt(B*B-A*d2*z3*z3));        % num. error possible - ok!
            if(mathjs.isComplex(z2) || mathjs.isNaN(z2) || !mathjs.isFinite(z2) || mathjs.smaller(z2, 0)) { //if ~isreal(z2) || isnan(z2) || isinf(z2) || z2 < 0 % num prob or wrong sign?
                if(mathjs.smaller(limit, -0.5)) { //if limit < -0.5                               % if we have no upper limit
                    z2 = mathjs.multiply(z1, mathjs.subtract(CONSTANT.EXT, 1)); //z2 = z1 * (EXT-1);                 % the extrapolate the maximum amount
                } else { //else
                    z2 = mathjs.divide(mathjs.subtract(limit, z1), 2);//z2 = (limit-z1)/2;                                   % otherwise bisect
                } //end
            } else if(mathjs.larger(limit, -0.5) && mathjs.complexLarger(mathjs.add(z2, z1), limit)) { //elseif (limit > -0.5) && (z2+z1 > limit)         % extraplation beyond max?
                z2 = mathjs.divide(mathjs.subtract(limit, z1), 2); //z2 = (limit-z1)/2;                                               % bisect
            } else if(mathjs.smaller(limit, -0.5) && mathjs.complexLarger(mathjs.add(z2, z1), mathjs.multiply(z1, CONSTANT.EXT))) { //elseif (limit < -0.5) && (z2+z1 > z1*EXT)       % extrapolation beyond limit
                z2 = mathjs.multiply(z1, mathjs.subtract(CONSTANT.EXT, 1.0)); //z2 = z1*(EXT-1.0);                           % set to extrapolation limit
            } else if(mathjs.smaller(z2, mathjs.multiply(-1, z3, CONSTANT.INT))) { //elseif z2 < -z3*INT
                z2 = mathjs.multiply(-1, z3, CONSTANT.INT); //z2 = -z3*INT;
            } else if(mathjs.larger(limit, -0.5) && mathjs.smaller(z2, mathjs.multiply(mathjs.subtract(limit, z1), mathjs.subtract(1.0, CONSTANT.INT)))) { //elseif (limit > -0.5) && (z2 < (limit-z1)*(1.0-INT))  % too close to limit?
                z2 = mathjs.multiply(mathjs.subtract(limit, z1), mathjs.subtract(1.0, CONSTANT.INT)); //z2 = (limit-z1)*(1.0-INT);
            } //end
            f3 = mathjs.clone(f2); //f3 = f2;
            d3 = mathjs.clone(d2); // d3 = d2;
            z3 = mathjs.multiply(-1, z2); // z3 = -z2;                  % set point 3 equal to point 2
            z1 = mathjs.add(z1, z2); //z1 = z1 + z2;
            z2ProductS = mathjs.multiply(z2, s);
            z2ProductS = mathjs.reshape(z2ProductS, [ z2ProductS.length, 1 ]);
            X = mathjs.add(X, z2ProductS); // X = X + z2*s;                      % update current estimates
            [ f2, df2 ] = updatedF(X); //[f2 df2] = eval(argstr);
            M = M - 1; //M = M - 1;
            i = i + ((length < 0)? 1 : 0); // i = i + (length<0);                             % count epochs?!
            d2 = mathjs.multiply(mathjs.transpose(df2), s); //d2 = df2'*s;
        } //end                                                      % end of line search
        var tmp;
        if(success) { //if success                                         % if line search succeeded
            f1 = mathjs.clone(f2); //f1 = f2;
            fX = mathjs.transpose(mathjs.concat(mathjs.transpose(fX), [ f1 ], 0)); // fX = [fX' f1]';
            console.log("SUCCESS", S, i, " | Cost: ", f1); //fprintf('%s %4i | Cost: %4.6e\r', S, i, f1);
            s = mathjs.subtract(
                    mathjs.multiply(
                        mathjs.dotDivide(
                            mathjs.subtract(
                                mathjs.multiply(mathjs.transpose(df2), df2),
                                mathjs.multiply(mathjs.transpose(df1), df2)
                            ),
                            mathjs.multiply(mathjs.transpose(df1), df1)
                        ),
                        s
                    ),
                    df2
                );
            //s = (df2'*df2-df1'*df2)/(df1'*df1)*s - df2;      % Polack-Ribiere direction
            tmp = mathjs.clone(df1); //tmp = df1;
            df1 = mathjs.clone(df2); // df1 = df2;
            df2 = mathjs.clone(tmp); // df2 = tmp;                         % swap derivatives
            d2 = mathjs.multiply(mathjs.transpose(df1), s); //d2 = df1'*s;
            if(mathjs.larger(d2, 0)) { //if d2 > 0                                      % new slope must be negative
                s = mathjs.multiply(-1, df1); //s = -df1;                              % otherwise use steepest direction
                d2 = mathjs.multiply(mathjs.multiply(-1, mathjs.transpose(s)), s);//d2 = -s'*s;
            } //end
            z1 = mathjs.multiply(z1, mathjs.min(CONSTANT.RATIO, mathjs.divide(d1, mathjs.subtract(d2, Number.MIN_VALUE)))); //z1 = z1 * min(RATIO, d1/(d2-realmin));          % slope ratio but max RATIO
            d1 = mathjs.clone(d2); //d1 = d2;
            lsFailed = 0; //ls_failed = 0;                              % this line search did not fail
        } else { //else
            console.log("FAILURE");
            X = mathjs.clone(X0); //X = X0;
            f1 = mathjs.clone(f0); // f1 = f0;
            df1 = mathjs.clone(df0); // df1 = df0;  % restore point from before failed line search
            if(lsFailed || mathjs.larger(i, mathjs.abs(length))) { //if ls_failed || i > abs(length)          % line search failed twice in a row
                break; //break;                             % or we ran out of time, so we give up
            } //end
            tmp = mathjs.clone(df1); //tmp = df1;
            df1 = mathjs.clone(df2); // df1 = df2;
            df2 = mathjs.clone(tmp); // df2 = tmp;                         % swap derivatives
            s = mathjs.multiply(-1, df1); //s = -df1;                                                    % try steepest
            d1 = mathjs.multiply(mathjs.multiply(-1, mathjs.transpose(s)), s); //d1 = -s'*s;
            z1 = mathjs.divide(1, mathjs.subtract(1, d1)); //z1 = 1/(1-d1);
            lsFailed = 1; //ls_failed = 1;                                    % this line search failed
        } //end
    } //end
    console.log("FMIN Conjugate Gradient Complete");
    return [ X, fX, i ];
};
