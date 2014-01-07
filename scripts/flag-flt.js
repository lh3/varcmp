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

var c, f_incl = 0, f_excl = 0;
while ((c = getopt(arguments, 'f:F:')) != null) {
	if (c == 'f') f_incl = parseInt(getopt.arg);
	else if (c == 'F') f_excl = parseInt(getopt.arg);
}

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t", 7);
	var m = /(0x[0-9a-fA-F]+)/.exec(t[6]);
	if (m == null) throw Error("No HEX filters");
	var f = parseInt(m[1]);
	if (f_incl && !(f&f_incl)) continue;
	if (f&f_excl) continue;
	print(buf);
}

buf.destroy();
file.close();
