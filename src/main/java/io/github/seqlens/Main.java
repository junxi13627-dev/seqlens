package io.github.seqlens;

import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * Command-line entry point.
 *
 * <pre>
 * java -jar seqlens.jar &lt;command&gt; [options] &lt;input.fasta | -&gt;
 * java -jar seqlens.jar &lt;command&gt; [options] --seq ATGC...
 * </pre>
 */
public final class Main {

    private static final String USAGE = """
            SeqLens 1.0.0 - DNA & protein sequence analysis

            Usage:
              seqlens <command> [options] <input.fasta>     read sequences from a FASTA file
              seqlens <command> [options] -                 read FASTA from standard input
              seqlens <command> [options] --seq <SEQUENCE>  analyze a raw sequence

            Commands:
              stats       length, base counts, GC content and melting temperature
              revcomp     reverse complement (FASTA output)
              transcribe  DNA -> mRNA (FASTA output)
              translate   DNA -> protein  [--frame 1|2|3, default 1]
              orfs        find open reading frames in all 6 frames  [--min-aa N, default 30]
              kmers       count k-mers  [--k N, default 3] [--top N, default 10]

            Examples:
              seqlens stats examples/sample.fasta
              seqlens orfs --min-aa 50 examples/sample.fasta
              seqlens translate --seq ATGGCCATTGTAATGGGCCGCTGA
            """;

    private Main() {
    }

    public static void main(String[] args) {
        int code = run(args);
        if (code != 0) {
            System.exit(code);
        }
    }

    /** Runs the CLI and returns an exit code (kept separate from main for testability). */
    static int run(String[] args) {
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            System.out.print(USAGE);
            return args.length == 0 ? 1 : 0;
        }
        try {
            Options opts = Options.parse(Arrays.copyOfRange(args, 1, args.length));
            List<FastaRecord> records = loadInput(opts);
            if (records.isEmpty()) {
                System.err.println("No sequences found in input.");
                return 1;
            }
            switch (args[0]) {
                case "stats" -> stats(records);
                case "revcomp" -> records.forEach(r ->
                        printFasta(r.id() + "_revcomp", DnaTools.reverseComplement(r.sequence())));
                case "transcribe" -> records.forEach(r ->
                        printFasta(r.id() + "_mRNA", DnaTools.transcribe(r.sequence())));
                case "translate" -> records.forEach(r ->
                        printFasta(r.id() + "_frame" + opts.frame,
                                Translator.translate(r.sequence(), opts.frame - 1)));
                case "orfs" -> orfs(records, opts.minAa);
                case "kmers" -> kmers(records, opts.k, opts.top);
                default -> {
                    System.err.println("Unknown command: " + args[0] + "\n");
                    System.err.print(USAGE);
                    return 2;
                }
            }
            return 0;
        } catch (IllegalArgumentException | IOException e) {
            System.err.println("Error: " + e.getMessage());
            return 2;
        }
    }

    private static List<FastaRecord> loadInput(Options opts) throws IOException {
        if (opts.rawSeq != null) {
            return List.of(new FastaRecord("input", "", opts.rawSeq));
        }
        if (opts.input == null) {
            throw new IllegalArgumentException("no input given (use a FASTA file, '-' or --seq)");
        }
        if (opts.input.equals("-")) {
            return FastaReader.read(new InputStreamReader(System.in, StandardCharsets.UTF_8));
        }
        return FastaReader.read(Path.of(opts.input));
    }

    private static void stats(List<FastaRecord> records) {
        System.out.printf("%-20s %10s %8s %8s %8s %8s %8s %8s %9s%n",
                "id", "length", "A", "C", "G", "T", "N", "GC%", "Tm(C)");
        for (FastaRecord r : records) {
            Map<Character, Integer> c = DnaTools.baseCounts(r.sequence());
            System.out.printf("%-20s %10d %8d %8d %8d %8d %8d %7.2f%% %9.1f%n",
                    truncate(r.id(), 20), r.length(),
                    c.get('A'), c.get('C'), c.get('G'), c.get('T'), c.get('N'),
                    DnaTools.gcContent(r.sequence()) * 100,
                    DnaTools.meltingTemperature(r.sequence()));
        }
    }

    private static void orfs(List<FastaRecord> records, int minAa) {
        System.out.printf("%-20s %6s %5s %8s %8s %8s  %s%n",
                "id", "strand", "frame", "start", "end", "aa_len", "protein");
        int total = 0;
        for (FastaRecord r : records) {
            for (Orf orf : OrfFinder.find(r.sequence(), minAa)) {
                System.out.printf("%-20s %6c %5d %8d %8d %8d  %s%n",
                        truncate(r.id(), 20), orf.strand(), orf.frame(), orf.start(), orf.end(),
                        orf.protein().length(), truncate(orf.protein(), 60));
                total++;
            }
        }
        System.out.println("\n" + total + " ORF(s) found with >= " + minAa + " amino acids.");
    }

    private static void kmers(List<FastaRecord> records, int k, int top) {
        for (FastaRecord r : records) {
            Map<String, Integer> counts = DnaTools.kmerCounts(r.sequence(), k);
            List<Map.Entry<String, Integer>> sorted = new ArrayList<>(counts.entrySet());
            sorted.sort(Map.Entry.<String, Integer>comparingByValue(Comparator.reverseOrder())
                    .thenComparing(Map.Entry.comparingByKey()));
            System.out.println("# " + r.id() + "  (k=" + k + ", " + counts.size() + " distinct)");
            sorted.stream().limit(top)
                    .forEach(e -> System.out.printf("%s\t%d%n", e.getKey(), e.getValue()));
        }
    }

    private static void printFasta(String header, String seq) {
        System.out.println(">" + header);
        for (int i = 0; i < seq.length(); i += 60) {
            System.out.println(seq.substring(i, Math.min(seq.length(), i + 60)));
        }
    }

    private static String truncate(String s, int max) {
        return s.length() <= max ? s : s.substring(0, max - 3) + "...";
    }

    /** Parsed command-line options. */
    static final class Options {
        String input;
        String rawSeq;
        int frame = 1;
        int minAa = 30;
        int k = 3;
        int top = 10;

        static Options parse(String[] args) {
            Options o = new Options();
            for (int i = 0; i < args.length; i++) {
                String a = args[i];
                switch (a) {
                    case "--seq" -> o.rawSeq = value(args, ++i, a).replaceAll("\\s", "");
                    case "--frame" -> o.frame = intValue(args, ++i, a, 1, 3);
                    case "--min-aa" -> o.minAa = intValue(args, ++i, a, 0, Integer.MAX_VALUE);
                    case "--k" -> o.k = intValue(args, ++i, a, 1, 32);
                    case "--top" -> o.top = intValue(args, ++i, a, 1, Integer.MAX_VALUE);
                    default -> {
                        if (a.startsWith("--")) {
                            throw new IllegalArgumentException("unknown option " + a);
                        }
                        o.input = a;
                    }
                }
            }
            return o;
        }

        private static String value(String[] args, int i, String flag) {
            if (i >= args.length) {
                throw new IllegalArgumentException(flag + " needs a value");
            }
            return args[i];
        }

        private static int intValue(String[] args, int i, String flag, int min, int max) {
            String v = value(args, i, flag);
            try {
                int n = Integer.parseInt(v);
                if (n < min || n > max) {
                    throw new IllegalArgumentException(flag + " must be between " + min + " and " + max);
                }
                return n;
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException(flag + " expects a number, got '" + v + "'");
            }
        }
    }
}
