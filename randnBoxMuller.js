module.exports = function randn_bm() {
    var u = 1 - Math.random();
    var v = 1 - Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
};
