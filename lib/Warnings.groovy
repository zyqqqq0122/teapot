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

    static String tolerance(Path tsv, String where) {
        def r = Files.readAllLines(tsv, StandardCharsets.ISO_8859_1)[1].split('\t', -1)
        "${where}: a light peptide and its heavy standard form a pair when measured within ${r[1]} s (${r[2]})".toString()
    }

    static String fragmentTolerance(Map t, String where, Object searchArgs) {
        def unit = Args.fragmentTolerance(searchArgs)[1] == 'th' ? 'Th' : 'ppm'
        (t.value != null
            ? "${where}: searched at a fragment tolerance of ${t.value} ${unit}, measured on the calibrants of ${t.runs.size()} runs together; a run without calibrants keeps -ftol"
            : "${where}: searched at -ftol of the search arguments (no run measured a tolerance, or search_tolerance off)").toString()
    }

    static String minCoelution(Path tsv, String where) {
        def r = Files.readAllLines(tsv, StandardCharsets.ISO_8859_1)[1].split('\t', -1)
        ((r[1] as double) > 0
            ? "${where}: a light peptide and its heavy standard form a pair when the light's signal follows the standard's over its peak (co-elution at least ${r[1]}, ${r[2]})"
            : "${where}: a light peptide and its heavy standard form a pair without a co-elution check (pair_min_coelution 0)").toString()
    }
}
