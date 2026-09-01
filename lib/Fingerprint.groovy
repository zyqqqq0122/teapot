import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.MessageDigest

class Fingerprint {

    static String of(p) {
        if (!p) return 'none'
        try {
            Path root = Paths.get(p.toString())
            if (!Files.exists(root)) return 'missing'
            def parts = []
            if (Files.isDirectory(root)) {
                Files.walk(root).withCloseable { s ->
                    s.filter { Files.isRegularFile(it) }.forEach { f ->
                        parts << "${root.relativize(f)}:${Files.size(f)}:${Files.getLastModifiedTime(f).toMillis()}".toString()
                    }
                }
            } else {
                parts << "${root.fileName}:${Files.size(root)}:${Files.getLastModifiedTime(root).toMillis()}".toString()
            }
            if (!parts) return 'empty'
            def md = MessageDigest.getInstance('MD5')
            md.update(parts.sort().join('|').getBytes('UTF-8'))
            return md.digest().encodeHex().toString()[0..11]
        }
        catch (Exception e) {
            return 'unreadable'
        }
    }
}
