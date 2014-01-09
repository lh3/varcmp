/******************
 *** From k8.js ***
 ******************/

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

/**********************
 *** Core algorithm ***
 **********************/

function cnv_cbs1(a0, start, end, winsize) // find the best cut for a single segment
{
	var a = a0.slice(start, end);
	a.sort(function(x,y){return x[1]-y[1]}); // sort by value
	var same = 1, T = 0;
	for (var i = 1; i <= a.length; ++i) { // conpute the rank
		if (i == a.length || a[i-1][1] != a[i][1]) {
			var rank = i - .5 * (same - 1);
			for (var j = i - same; j <= i - 1; ++j)
				a[j][2] = rank;
			if (same > 1) T += (same - 1.) * same * (same + 1.), same = 1; // T to correct for ties
		} else ++same;
	}
	var z = ((a.length + 1) - T / a.length / (a.length - 1)) / 12.;
	a.sort(function(x,y){return x[0]-y[0]});
	var s = [], sum = 0;
	s[0] = 0;
	for (var i = 0; i < a.length; ++i) s[i+1] = s[i] + a[i][2];
	var max = [0, -1, -1, 0];
	for (var j = 1; j <= a.length; ++j) {
		for (var i = j > winsize? j - winsize : 0; i < j; ++i) { // compute the standardized Mann-Whitney U
			var n = j - i, m = a.length - n;
			var x = (s[j] - s[i] - .5 * n * (a.length + 1)) / Math.sqrt(m * n * z);
			var y = x > 0? x : -x;
			if (y > max[0]) max = [y, i, j, x]; // this is the best cut
		}
	}
	return max;
}

function cnv_cbs(chr, a, thres, winsize) // recursively cut
{ // FIXME: it does not consider the edge effect; should not be a big issue
	function print_intv(b, e) { // print an interval
		if (b == e) return;
		var t = a.slice(b, e);
		t.sort(function(x,y){return x[1]-y[1]}); // sort to compute the median
		print(chr, t[0][0], t[t.length-1][0] + 1, t.length, t[Math.floor(t.length/2)][1]);
	}

	warn("Processing chromsome " + chr + "...");
	var queue = [[0, a.length, -1]];
	while (queue.length) {
		var e = queue.shift();
		var start = e[0], end = e[1];
		var cut = cnv_cbs1(a, start, end, winsize);
		if (cut[0] < -thres || cut[0] > thres) { // a significant cut; then cut further
			if (cut[1] >= 4) queue.push([start, start + cut[1]]);
			else print_intv(start, start + cut[1]);
			if (cut[2] - cut[1] >= 4) queue.push([start + cut[1], start + cut[2]]);
			else print_intv(start + cut[1], start + cut[2]);
			if (end - start - cut[2] >= 4) queue.push([start + cut[2], end]);
			else print_intv(start + cut[2], end);
		} else print_intv(start, end); // no significant cut can be found
	}
}

/*********************
 *** Main function ***
 *********************/

function main(args)
{
	var c, thres = 3., winsize = 100;
	while ((c = getopt(args, 't:w:')) != null) {
		if (c == 't') thres = parseFloat(getopt.arg);
		else if (c == 'w') winsize = parseInt(getopt.arg);
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var last_chr = null, a = [];

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] != last_chr) {
			if (a.length) cnv_cbs(last_chr, a, thres, winsize);
			last_chr = t[0];
			a = [];
		}
		a.push([parseInt(t[1]), parseFloat(t[2])]);
	}
	if (a.length) cnv_cbs(last_chr, a, thres, winsize);

	buf.destroy();
	file.close();
}

main(arguments);
