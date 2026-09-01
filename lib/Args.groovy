class Args {
	static String flat(Object value) {
		return value == null ? '' : value.toString().replaceAll(/\s+/, ' ').trim()
	}

	static final List<String> INHERITED_FLAGS = [
		'-ftol', '-ftolunits',
		'-lftol', '-lftolunits',
		'-ptol', '-ptolunits',
		'-foffset', '-poffset',
		'-frag',
		'-fixed',
		'-enzyme',
	]

	static String inherit(Object searchArgs) {
		def toks = flat(searchArgs).split(' ').findAll { it }
		def out = []
		for (int i = 0; i < toks.size(); i++) {
			if (INHERITED_FLAGS.contains(toks[i]) && i + 1 < toks.size() && !toks[i + 1].startsWith('-')) {
				out << toks[i] << toks[i + 1]
				i++
			}
		}
		return out.join(' ')
	}
}
