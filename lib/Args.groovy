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

	static String mode(Object value) {
		if (value in [false, 'false', 'off', 'OFF', 'Off']) return 'off'
		if (value in [true, 'true', 'auto']) return 'auto'
		return null
	}

	static List fragmentTolerance(Object searchArgs) {
		def toks = flat(searchArgs).split(' ').findAll { it }
		def opt = [:]
		for (int i = 0; i + 1 < toks.size(); i++) {
			if (toks[i].startsWith('-')) opt[toks[i]] = toks[i + 1]
		}
		def tol = (opt['-ftol'] ?: '10') as double
		def unit = (opt['-ftolunits'] ?: 'ppm').toLowerCase()
		return unit == 'amu' ? [tol, 'th'] : unit == 'resolution' ? [0.5e6 / tol, 'ppm'] : [tol, 'ppm']
	}

	static String withFragmentTolerance(Object searchArgs, Object tol) {
		def unit = fragmentTolerance(searchArgs)[1] == 'th' ? 'AMU' : 'PPM'
		def toks = flat(searchArgs).split(' ').findAll { it }
		def put = { String flag, String v ->
			def i = toks.indexOf(flag)
			if (i >= 0 && i + 1 < toks.size()) toks[i + 1] = v
			else toks += [flag, v]
		}
		put('-ftol', tol.toString())
		if (toks.contains('-ftolunits') || unit != 'PPM') put('-ftolunits', unit)
		return toks.join(' ')
	}

	static String withBatchTolerance(Object searchArgs, Map t) {
		return t.value != null ? withFragmentTolerance(searchArgs, t.value) : flat(searchArgs)
	}

	static Map readTolerance(Object tsv, Object runs) {
		def lines = java.nio.file.Files.readAllLines(tsv as java.nio.file.Path)
		def head = lines[0].split('\t', -1), cells = lines[1].split('\t', -1)
		def at = cells[head.findIndexOf { it == 'runs_at_it' }] as int
		def rows = java.nio.file.Files.readAllLines(runs as java.nio.file.Path)
		def cols = rows[0].split('\t', -1) as List
		def own = rows.drop(1).findAll { it }.collect { it.split('\t', -1) }
			.findAll { it[cols.indexOf('source')] == 'measured' }
			.collectEntries { [(it[0]): it[cols.indexOf('searched_at')] as double] }
		return [value: at == 0 ? null : cells[head.findIndexOf { it == 'search_tolerance' }] as double, runs: own]
	}

	static String withRunTolerance(Object searchArgs, Map t, String sid) {
		return t.runs.containsKey(sid) ? withFragmentTolerance(searchArgs, t.runs[sid]) : flat(searchArgs)
	}

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
