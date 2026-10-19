import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.context.IsolationWindow;
import edu.washington.gs.maccoss.encyclopedia.context.IsolationWindowReader;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class MassListDecoyGenerator {

    public static void main(String[] args) throws Exception {
        String in = null, out = null;
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals("-massList")) in = args[i + 1];
            if (args[i].equals("-o"))        out = args[i + 1];
        }
        if (in == null || out == null) {
            System.err.println("usage: MassListDecoyGenerator -massList <file> -o <file>");
            System.exit(1);
        }
        File input = new File(in);
        if (!input.exists()) {
            System.err.println("Mass list not found: " + input.getAbsolutePath());
            System.exit(1);
        }
        ArrayList<IsolationWindow> windows = addDecoys(input.getAbsolutePath());
        write(windows, new File(out));
    }

    public static boolean hasDecoys(String massListPath) {
        for (IsolationWindow w : IsolationWindowReader.parseMassList(massListPath)) {
            if (w.isDecoy()) return true;
        }
        return false;
    }

    public static ArrayList<IsolationWindow> addDecoys(String massListPath) {
        ArrayList<IsolationWindow> input = IsolationWindowReader.parseMassList(massListPath);

        for (IsolationWindow w : input) {
            if (w.isDecoy()) {
                System.err.println("Mass list already contains decoys; leaving it unchanged");
                return input;
            }
        }

        AminoAcidConstants constants = new AminoAcidConstants();

        HashSet<String> takenSequences = new HashSet<String>();
        ArrayList<IsolationWindow> output = new ArrayList<IsolationWindow>();

        for (IsolationWindow target : input) {
            output.add(target);

            String sequence = target.getCompound();
            byte charge = target.getCharge();
            takenSequences.add(sequence);

            String decoy = reverseText(sequence);
            String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
            double decoyMz = constants.getChargedMass(correctedDecoyMass, charge);
            takenSequences.add(decoy);

            output.add(new IsolationWindow(decoy, decoyMz, charge,
                                           target.getRtMin(), target.getRtMax(), true));
        }
        return output;
    }

    static String reverseText(String sequence) {
        ArrayList<String> aas = new ArrayList<String>();
        String nterm = "";
        int i = 0;
        while (i < sequence.length()) {
            char c = sequence.charAt(i);
            if (c == '[') {
                int end = sequence.indexOf(']', i);
                String mod = sequence.substring(i, end + 1);
                if (aas.isEmpty()) {
                    nterm += mod;
                } else {
                    aas.set(aas.size() - 1, aas.get(aas.size() - 1) + mod);
                }
                i = end + 1;
            } else {
                aas.add(aas.isEmpty() ? c + nterm : Character.toString(c));
                i++;
            }
        }
        for (int start = 1, stop = aas.size() - 2; start < stop; start++, stop--) {
            String t = aas.get(start);
            aas.set(start, aas.get(stop));
            aas.set(stop, t);
        }
        return String.join("", aas);
    }

    private static void write(ArrayList<IsolationWindow> windows, File output) throws IOException {
        Path parent = output.getAbsoluteFile().toPath().getParent();
        if (parent != null) Files.createDirectories(parent);

        int decoys = 0;
        try (BufferedWriter writer = Files.newBufferedWriter(output.toPath())) {
            writer.write("Compound\tFormula\tAdduct\tm/z\tz\tRT Time (min)\tWindow (min)\tisDecoy");
            writer.newLine();
            for (IsolationWindow window : windows) {
                if (window.isDecoy()) decoys++;
                float rtCenterMin = ((window.getRtMin() + window.getRtMax()) / 2.0f) / 60.0f;
                float windowMin   = (window.getRtMax() - window.getRtMin()) / 60.0f;
                writer.write(window.getCompound() + "\t" + "\t" + "(no adduct)" + "\t"
                        + window.getTargetMz() + "\t" + window.getCharge() + "\t"
                        + rtCenterMin + "\t" + windowMin + "\t" + window.isDecoy());
                writer.newLine();
            }
        }
        System.err.println("Wrote " + output.getAbsolutePath() + ": "
                + (windows.size() - decoys) + " targets, " + decoys + " decoys");
    }
}
