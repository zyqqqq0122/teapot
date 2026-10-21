import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path

class Warnings {

    static List<String> inFile(Path file, String where) {
        Files.readAllLines(file, StandardCharsets.ISO_8859_1).findAll { it.trim() }
                                .collect { "${where}: ${it.trim()}".toString() }
    }

    static List<String> inLog(Path log, String where) {
        Files.readAllLines(log, StandardCharsets.ISO_8859_1).findAll { it.contains('WARNING') }
                               .collect { "${where}: ${it.replaceFirst(/^.*?WARNING:?\s*/, '').trim()}".toString() }
    }

    static List<Double> perRun(Path tsv, String column) {
        def lines = Files.readAllLines(tsv, StandardCharsets.ISO_8859_1)
        def i = lines[0].split('\t', -1).findIndexOf { it == column }
        lines.drop(1).findAll { it.trim() && !it.startsWith('*\t') }.collect { it.split('\t', -1)[i] as Double }
    }

    static String span(List<Double> v, String unit) {
        v.min() == v.max() ? "${v.min()}${unit}" : "${v.min()} to ${v.max()}${unit} (each run's own)"
    }

    static String tolerance(Path tsv, String where) {
        def v = perRun(tsv, 'tolerance_seconds')
        (v ? "${where}: a light peptide and its heavy standard form a pair when measured within ${span(v, ' s')}, half the median peak width of the standards of the run and its neighbours"
           : "${where}: no run to derive a pair tolerance from").toString()
    }

    static String fragmentTolerance(Map t, String where, Object searchArgs) {
        def unit = Args.fragmentTolerance(searchArgs)[1] == 'th' ? 'Th' : 'ppm'
        def v = t.runs.values() as List
        (v ? "${where}: searched at a fragment tolerance of ${span(v, ' ' + unit)}, measured on the calibrants of each run and its neighbours; a run without calibrants keeps -ftol"
           : "${where}: searched at -ftol of the search arguments (no run measured a tolerance, or search_tolerance off)").toString()
    }

    static String minCoelution(Path tsv, String where) {
        def v = perRun(tsv, 'min_coelution')
        (!v ? "${where}: no run measured at the standard's peak"
            : v.max() > 0 ? "${where}: a light peptide and its heavy standard form a pair when the light's signal follows the standard's over its peak (co-elution at least ${span(v, '')}; qc in quant/${where}/pair_min_coelution.tsv)"
            : "${where}: a light peptide and its heavy standard form a pair without a co-elution check (pair_min_coelution 0)").toString()
    }
}
