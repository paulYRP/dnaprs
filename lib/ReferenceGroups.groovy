import java.nio.file.Files
import java.nio.file.FileVisitOption
import java.nio.file.Path

class ReferenceGroups {
    static List balance(List<Path> roots) {
        def paths = []
        roots.each { root ->
            if (Files.isDirectory(root)) {
                Files.walk(root, 2, FileVisitOption.FOLLOW_LINKS).withCloseable { stream ->
                    stream.filter { Files.isRegularFile(it) }.forEach { paths << it }
                }
            } else {
                paths << root
            }
        }
        def records = paths.unique().findResults { source ->
            def match = source.fileName.toString() =~ /^chr([1-9]|1[0-9]|2[0-2])(?:[._].*)?\.(?:bref3|vcf(?:\.gz|\.bgz)?|bcf)$/
            match.matches() ? [chromosome: match[0][1] as Integer, path: source, bytes: Files.size(source)] : null
        }
        if (!records) throw new IllegalArgumentException('Reference panel contains no supported autosomal chromosome files.')
        def repeated = records.groupBy { it.chromosome }.findAll { chromosome, values -> values.size() > 1 }.keySet().sort()
        if (repeated) throw new IllegalArgumentException("Reference panel contains repeated chromosome files: ${repeated.join(', ')}.")
        def groups = (1..Math.min(6, records.size())).collect { [id: it, bytes: 0L, records: []] }
        records.sort { a, b -> (b.bytes <=> a.bytes) ?: (a.chromosome <=> b.chromosome) }.each { record ->
            def group = groups.min { a, b -> (a.bytes <=> b.bytes) ?: (a.id <=> b.id) }
            group.records << record
            group.bytes += record.bytes
        }
        groups.each { it.records.sort { record -> record.chromosome } }
        groups
    }
}
