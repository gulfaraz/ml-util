module.exports = function randomPermutation(n) {
    var i;
    var numbers = [];
    for (i = 0; i < n; i++) {
        numbers[i] = i;
    }
    for (i = 0; i < n; i++) {
        var pos = n - i - 1;
        var spos = parseInt(Math.random() * (pos + 1));
        var tmp = numbers[spos];
        numbers[spos] = numbers[pos];
        numbers[pos] = tmp;
    }
    return numbers
};
