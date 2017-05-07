var mathjs = require("mathjs");

mathjs.isComplex = function isComplex(value) {
    return value != null && value.isComplex;
};

mathjs.complexLarger = function complexLarger(a, b) {
    var isLarger = true;
    if(mathjs.isComplex(a) || mathjs.isComplex(b)) {
        isLarger = mathjs.abs(a) > mathjs.abs(b)
        || (mathjs.abs(a) == mathjs.abs(b)
            && mathjs.arg(b) > mathjs.arg(b));
    } else {
        isLarger = mathjs.larger(a, b);
    }
    return isLarger;
};

mathjs.complexSmaller = function complexSmaller(a, b) {
    var isSmaller = true;
    if(mathjs.isComplex(a) || mathjs.isComplex(b)) {
        isSmaller = mathjs.abs(a) < mathjs.abs(b)
        || (mathjs.abs(a) == mathjs.abs(b)
            && mathjs.arg(b) < mathjs.arg(b));
    } else {
        isSmaller = mathjs.smaller(a, b);
    }
    return isSmaller;
};

mathjs.isFinite = function isFinite(value) {
  if (value == Number.POSITIVE_INFINITY
    || value == Number.NEGATIVE_INFINITY){
    return false;
  }
  if (value instanceof Array){
    for (var i=0; i<value.length; i++){
      var item=value[i];
      if (item == Number.POSITIVE_INFINITY
        || item == Number.NEGATIVE_INFINITY
        || math.equals(item, Number.POSITIVE_INFINITY)
        || math.equals(item, Number.NEGATIVE_INFINITY)){
        return false;
      }
      return true;
    }
  }
  return true;
};

module.exports = function fmincg(f, X, options, P1, P2, P3, P4, P5) {
    console.log("Starting FMIN Conjugate Gradient Function...");

    var length = 100;

    if(options && !isNaN(options.maxIterations)) {
        length = options.maxIterations;
    }

    var CONSTANT = {
        "RHO" : 0.01,
        "SIG" : 0.5,
        "INT" : 0.1,
        "EXT" : 3.0,
        "MAX" : 20,
        "RATIO" : 100
    };

    var argumentList = Array.prototype.slice.call(arguments);
    argumentList.splice(0, 3);

    var updatedF = function updatedF(X) {
        var updatedArgumentList = argumentList.slice(0);
        updatedArgumentList.unshift(X);
        return f.apply(null, updatedArgumentList);
    };

    var red = 1;
    var S = "Iteration - ";

    if(length && length.length === 2) {
        length = length[1];
        red = length[2];
    }

    var i = 0;
    var lsFailed = 0;
    var fX = [];

    var [ f1, df1 ] = updatedF(X);

    i = i + ((length<0) ? 1 : 0);

    var s = mathjs.multiply(-1, df1);
    var d1 = mathjs.multiply(
        mathjs.multiply(-1, mathjs.transpose(s)),
        s
    );
    var z1 = mathjs.dotDivide(red, mathjs.subtract(1, d1));

    while(i < mathjs.abs(length)) {
        i = i + ((length>0) ? 1 : 0);
        var X0 = mathjs.clone(X);
        var f0 = mathjs.clone(f1);
        var df0 = mathjs.clone(df1);
        var z1ProductS = mathjs.multiply(s, z1);
        var z1ProductSLength = 0;
        mathjs.forEach(
            z1ProductS,
            function countLengthOfZ1ProductS(value, index, matrix) {
                z1ProductSLength++;
            }
        );
        z1ProductS = mathjs.reshape(
            z1ProductS,
            [ z1ProductSLength, 1 ]
        );
        X = mathjs.add(X, z1ProductS);
        var [ f2, df2 ] = updatedF(X);
        i = i + ((length<0) ? 1 : 0);
        var d2 = mathjs.multiply(mathjs.transpose(df2), s);
        var f3 = mathjs.clone(f1);
        var d3 = mathjs.clone(d1);
        var z3 = mathjs.multiply(-1, z1);
        var M = mathjs.min(
            CONSTANT.MAX,
            mathjs.subtract(mathjs.multiply(-1, length), i)
        );
        if(length > 0) {
            M = CONSTANT.MAX;
        }
        var success = 0;
        var limit = -1;
        while(1) {
            while(
                (
                    mathjs.larger(
                        f2,
                        mathjs.add(
                            f1,
                            mathjs.multiply(
                                mathjs.multiply(z1, CONSTANT.RHO),
                                d1
                            )
                        )
                    )
                    || mathjs.larger(
                        d2,
                        mathjs.multiply(
                            mathjs.multiply(-1, CONSTANT.SIG),
                            d1
                        )
                    )
                )
                && (M > 0)
            ) {
                limit = mathjs.clone(z1);
                var z2;
                if(mathjs.larger(f2, f1)) {
                    z2 = mathjs.subtract(
                        z3,
                        mathjs.divide(
                            mathjs.multiply(
                                mathjs.multiply(
                                    mathjs.multiply(0.5, d3),
                                    z3
                                ),
                                z3
                            ),
                            mathjs.add(
                                mathjs.subtract(
                                    mathjs.multiply(d3, z3),
                                    f3
                                ),
                                f2
                            )
                        )
                    );
                } else {
                    var A = mathjs.add(
                        mathjs.multiply(
                            6,
                            mathjs.divide(
                                mathjs.subtract(f2,f3),
                                z3
                            )
                        ),
                        mathjs.multiply(
                            3,
                            mathjs.add(d2,d3)
                        )
                    );
                    var B = mathjs.subtract(
                        mathjs.multiply(
                            3,
                            mathjs.subtract(f3, f2)
                        ),
                        mathjs.multiply(
                            z3,
                            mathjs.add(
                                d3,
                                mathjs.multiply(2, d2)
                            )
                        )
                    );
                    z2 = mathjs.divide(
                        mathjs.subtract(
                            mathjs.sqrt(
                                mathjs.subtract(
                                    mathjs.multiply(B, B),
                                    mathjs.multiply(
                                        mathjs.multiply(
                                            mathjs.multiply(A, d2),
                                            z3
                                        ),
                                        z3
                                    )
                                )
                            ),
                            B
                        ),
                        A
                    );
                }
                if(mathjs.isNaN(z2) || !mathjs.isFinite(z2)) {
                    z2 = mathjs.divide(z3, 2);
                }
                z2 = mathjs.max(
                    mathjs.min(
                        mathjs.sum(z2),
                        mathjs.sum(
                            mathjs.multiply(CONSTANT.INT, z3)
                        )
                    ),
                    mathjs.sum(
                        mathjs.multiply(
                            mathjs.subtract(1, CONSTANT.INT),
                            z3
                        )
                    )
                );
                z1 = mathjs.add(z1, z2);
                var z2ProductS = mathjs.multiply(z2, s);
                var z2ProductSLength = 0;
                mathjs.forEach(
                    z2ProductS,
                    function countLengthOfZ2ProductS(
                        value,
                        index,
                        matrix
                    ) {
                        z2ProductSLength++;
                    }
                );
                z2ProductS = mathjs.reshape(
                    z2ProductS,
                    [ z2ProductSLength, 1 ]
                );
                X = mathjs.add(X, z2ProductS);
                [ f2, df2 ] = updatedF(X);
                M = M - 1;
                i = i + ((length<0) ? 1 : 0);
                d2 = mathjs.multiply(mathjs.transpose(df2), s);
                z3 = mathjs.subtract(z3, z2);
            }
            if(
                mathjs.larger(f2,
                    mathjs.add(f1,
                        mathjs.multiply(
                            mathjs.multiply(z1, CONSTANT.RHO),
                            d1
                        )
                    )
                )
                || mathjs.larger(
                    d2,
                    mathjs.multiply(
                        mathjs.multiply(-1, CONSTANT.SIG),
                        d1
                    )
                )
            ) {
                break;
            } else if(
                mathjs.larger(
                    d2,
                    mathjs.multiply(CONSTANT.SIG, d1)
                )
            ) {
                success = 1;
                break;
            } else if(M === 0) {
                break;
            }
            var z2;
            var A = mathjs.add(
                mathjs.multiply(
                    6,
                    mathjs.divide(
                        mathjs.subtract(f2, f3),
                        z3
                    )
                ),
                mathjs.multiply(
                    3,
                    mathjs.add(d2, d3)
                )
            );
            var B = mathjs.subtract(
                mathjs.multiply(
                    3,
                    mathjs.subtract(f3, f2)
                ),
                mathjs.multiply(
                    z3,
                    mathjs.add(
                        d3,
                        mathjs.multiply(2, d2)
                    )
                )
            );
            z2 = mathjs.divide(
                mathjs.multiply(
                    mathjs.multiply(
                        mathjs.multiply(-1, d2),
                        z3
                    ),
                    z3
                ),
                mathjs.add(
                    B,
                    mathjs.sqrt(
                        mathjs.subtract(
                            mathjs.multiply(B, B),
                            mathjs.multiply(
                                mathjs.multiply(
                                    mathjs.multiply(A, d2),
                                    z3
                                ),
                                z3
                            )
                        )
                    )
                )
            );
            if(
                mathjs.isComplex(z2)
                || mathjs.isNaN(z2)
                || !mathjs.isFinite(z2)
                || mathjs.complexSmaller(z2, 0)
            ) {
                if(mathjs.complexSmaller(limit, -0.5)) {
                    z2 = mathjs.multiply(
                        z1,
                        mathjs.subtract(CONSTANT.EXT, 1)
                    );
                } else {
                    z2 = mathjs.divide(
                        mathjs.subtract(limit, z1),
                        2
                    );
                }
            } else if(
                mathjs.larger(limit, -0.5)
                && mathjs.complexLarger(
                    mathjs.add(z2, z1),
                    limit
                )
            ) {
                z2 = mathjs.divide(
                    mathjs.subtract(limit, z1),
                    2
                );
            } else if(
                mathjs.complexSmaller(limit, -0.5)
                && mathjs.complexLarger(
                    mathjs.add(z2, z1),
                    mathjs.multiply(z1, CONSTANT.EXT)
                )
            ) {
                z2 = mathjs.multiply(
                    z1,
                    mathjs.subtract(CONSTANT.EXT, 1.0)
                );
            } else if(
                mathjs.complexSmaller(
                    z2,
                    mathjs.multiply(
                        mathjs.multiply(-1, z3),
                        CONSTANT.INT
                    )
                )
            ) {
                z2 = mathjs.multiply(
                    mathjs.multiply(-1, z3),
                    CONSTANT.INT
                );
            } else if(
                mathjs.larger(limit, -0.5)
                && mathjs.complexSmaller(
                    z2,
                    mathjs.multiply(
                        mathjs.subtract(limit, z1),
                        mathjs.subtract(1.0, CONSTANT.INT)
                    )
                )
            ) {
                z2 = mathjs.multiply(
                    mathjs.subtract(limit, z1),
                    mathjs.subtract(1.0, CONSTANT.INT)
                );
            }
            f3 = mathjs.clone(f2);
            d3 = mathjs.clone(d2);
            z3 = mathjs.multiply(-1, z2);
            z1 = mathjs.add(z1, z2);
            z2ProductS = mathjs.multiply(z2, s);
            z2ProductSLength = 0;
            mathjs.forEach(
                z2ProductS,
                function countLengthOfZ2ProductS(
                    value,
                    index,
                    matrix
                ) {
                    z2ProductSLength++;
                }
            );
            z2ProductS = mathjs.reshape(
                z2ProductS,
                [ z2ProductSLength, 1 ]
            );
            X = mathjs.add(X, z2ProductS);
            [ f2, df2 ] = updatedF(X);
            M = M - 1;
            i = i + ((length < 0)? 1 : 0);
            d2 = mathjs.multiply(mathjs.transpose(df2), s);
        }
        var tmp;
        if(success) {
            f1 = mathjs.clone(f2);
            fX = mathjs.transpose(
                mathjs.concat(mathjs.transpose(fX),[ f1 ], 0)
            );
            console.log("SUCCESS", S, i, " | Cost: ", f1);
            s = mathjs.subtract(
                    mathjs.multiply(
                        mathjs.dotDivide(
                            mathjs.subtract(
                                mathjs.multiply(
                                    mathjs.transpose(df2),
                                    df2
                                ),
                                mathjs.multiply(
                                    mathjs.transpose(df1),
                                    df2
                                )
                            ),
                            mathjs.multiply(
                                mathjs.transpose(df1),
                                df1
                            )
                        ),
                        s
                    ),
                    df2
                );
            tmp = mathjs.clone(df1);
            df1 = mathjs.clone(df2);
            df2 = mathjs.clone(tmp);
            d2 = mathjs.multiply(mathjs.transpose(df1), s);
            if(mathjs.larger(d2, 0)) {
                s = mathjs.multiply(-1, df1);
                d2 = mathjs.multiply(
                    mathjs.multiply(
                        -1,
                        mathjs.transpose(s)
                    ),
                    s
                );
            }
            z1 = mathjs.multiply(
                z1,
                mathjs.min(
                    CONSTANT.RATIO,
                    mathjs.divide(
                        d1,
                        mathjs.subtract(d2, Number.MIN_VALUE)
                    )
                )
            );
            d1 = mathjs.clone(d2);
            lsFailed = 0;
        } else {
            console.log("FAILURE");
            X = mathjs.clone(X0);
            f1 = mathjs.clone(f0);
            df1 = mathjs.clone(df0);
            if(lsFailed || mathjs.larger(i, mathjs.abs(length))) {
                break;
            }
            tmp = mathjs.clone(df1);
            df1 = mathjs.clone(df2);
            df2 = mathjs.clone(tmp);
            s = mathjs.multiply(-1, df1);
            d1 = mathjs.multiply(
                athjs.multiply(
                    -1,
                    mathjs.transpose(s)
                ),
                s
            );
            z1 = mathjs.divide(1, mathjs.subtract(1, d1));
            lsFailed = 1;
        }
    }
    console.log("FMIN Conjugate Gradient Complete");
    return [ X, fX, i ];
};
