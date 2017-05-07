var mathjs = require("mathjs");

module.exports = function sigmoid(z) {
    var g = mathjs.dotDivide(
        1.0,
        mathjs.add(
            1.0,
            mathjs.exp(
                mathjs.multiply(-1, z)
            )
        )
    );
    return g;
};
