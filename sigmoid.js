var mathjs = require("mathjs");

module.exports = function (z) {
    console.log("Starting Sigmoid Function...");

    var g = mathjs.dotDivide(1.0, mathjs.add(1.0, mathjs.exp(mathjs.multiply(-1, z))));

    console.log("Sigmoid Function Complete");
    return g;
};
