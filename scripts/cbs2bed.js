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

var c, bound = 3.;
while ((c = getopt(arguments, 'b:')) != null)
	if (c == 'b') bound = parseFloat(getopt.arg);

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();
var a = [];

while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	for (var i = 1; i < t.length; ++i)
		t[i] = parseInt(t[i]);
	t[5] = t[4] + Math.random() - .5;
	a.push(t);
}

a.sort(function(x,y){return x[5]-y[5]});

var b = [];
for (var i = 0; i < a.length; ++i)
	for (var j = 0; j < a[i][3]; ++j)
		b.push(a[i][5]);

var p25 = b[Math.floor(b.length * .25)];
var p75 = b[Math.floor(b.length * .75)];
var cutoff = p75 + bound * (p75 - p25);
warn(p25, p75, cutoff);

for (var i = 0; i < a.length; ++i)
	if (a[i][4] > cutoff) print(a[i].slice(0, 3).join("\t"));

buf.destroy();
file.close();
