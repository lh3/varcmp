var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

//var c, AB = 30, LC = 0, DP_coef = 3., DS = 1, FS = -20, min_q = 30, min_dp = 0;
var c, AB = 20, LC = 0, DP_coef = 4, DS = 0, FS = -30, min_q = 30, min_dp = 8;
while ((c = getopt(arguments, "La:l:d:D:F:q:")) != null) {
	if (c == 'a') AB = parseInt(getopt.arg);
	else if (c == 'q') min_q = parseFloat(getopt.arg);
	else if (c == 'L') LC = -1000000000;
	else if (c == 'F') FS = parseFloat(getopt.arg);
}

if (getopt.ind + 1 > arguments.length) {
	print("Usage: k8 vcf-extra-flt.js <var-extra.vcf>");
	exit(1);
}

var file = new File(arguments[getopt.ind]);
var buf = new Bytes();

// pass 1: compute avg depth
var sum_dp = 0, n_dp = 0;
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	var m;
	if (parseFloat(t[5]) < min_q) continue;
	if ((m = /qLC=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < LC) continue;
	if ((m = /qDP=(-?\d+)/.exec(t[7])) == null) continue;
	var dp = -parseInt(m[1]);
	sum_dp += dp; ++n_dp;
}
var avg_dp = sum_dp / n_dp;
var max_dp = Math.sqrt(avg_dp) * DP_coef + avg_dp;
var DP = -max_dp;
file.close();
warn(n_dp, avg_dp, max_dp);

// pass 2: filter
file = new File(arguments[getopt.ind]);
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	var m;
	if (parseFloat(t[5]) < min_q) continue;
	if ((m = /qLC=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < LC) continue;
	if ((m = /qDP=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < DP) continue;
	if ((m = /qDP=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) > -min_dp) continue;
	if ((m = /qDS=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < DS) continue;
	if ((m = /qFS=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < FS) continue;
	if ((m = /qAB=(-?\d+)/.exec(t[7])) != null && parseInt(m[1]) < AB) continue;
	print(buf);
}
file.close();

buf.destroy();
