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

    static String tolerance(Path tsv, String where, Object primary, Object minCoelution) {
        if (primary == 'standard_peak') {
            return ((minCoelution as double) > 0
                ? "${where}: a light peptide and its heavy standard form a pair when the light's signal follows the standard's over its peak (co-elution at least ${minCoelution})"
                : "${where}: a light peptide and its heavy standard form a pair without a co-elution check (pair_min_coelution 0)").toString()
        }
        def r = Files.readAllLines(tsv, StandardCharsets.ISO_8859_1)[1].split('\t', -1)
        "${where}: a light peptide and its heavy standard form a pair when measured within ${r[1]} s (${r[2]})".toString()
    }
}
