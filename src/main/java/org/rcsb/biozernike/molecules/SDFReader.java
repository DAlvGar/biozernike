package org.rcsb.biozernike.molecules;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

public class SDFReader {

    public static Point3d[] getAtomPositions(IAtomContainer molecule) {        
        Point3d[] points = new Point3d[molecule.getAtomCount()];
		for(int i = 0; i< molecule.getAtomCount();i++) {
			points[i] = molecule.getAtom(i).getPoint3d();
		}
        return points;
    }

    public static Map<String, Object> parseString(String inputString) {
        Map<String, Object> resultMap = new HashMap<>();
        String patternString = "(?<molid>[A-Za-z0-9]+)-(?<active>[A-Za-z]+[A-Za-z0-9_]*)\\.(?<conformer>\\d+)";
        Pattern pattern = Pattern.compile(patternString);
        Matcher matcher = pattern.matcher(inputString);

        if (matcher.matches()) {
            String molid = matcher.group("molid");
            String active_s = matcher.group("active");
            int conformer = Integer.parseInt(matcher.group("conformer"));
            Boolean active = false;
            // Special handling for 'active' field
            if (active_s.contains("_")) {
                active_s = active_s.split("_")[0];
                conformer += 100; // avoid repeated conformer numbers
            }
            if (active_s.contains("dec")) {
               active = false;
            } else {
               active = true;
            }

            // Populate the result map
            resultMap.put("molid", molid);
            resultMap.put("active", active);
            resultMap.put("conformer", conformer);
            resultMap.put("original", inputString);
        } else {
            System.out.println("ERROR WITH: " + inputString);
        }

        return resultMap;
    }

    public static List<List<Double>> getAtomFields(IAtomContainer molecule) {
        List<List<Double>> fields = new ArrayList<>();

        // Define a regular expression to match floating-point numbers
        String regex = "[-+]?\\d*\\.?\\d+";
        // Create a Pattern object
        Pattern pattern = Pattern.compile(regex);
        for (IAtom atom : molecule.atoms()) {
            List<Double> afields = new ArrayList<Double>();
            Map<Object, Object> p = atom.getProperties();
            String entry = p.entrySet().toString();
            // Create a Matcher object
            Matcher matcher = pattern.matcher(entry);
            // List to store extracted numbers
            int row_count=0;
            while (matcher.find() & row_count < 5) {
                // Convert the matched string to Double and add to the list
                afields.add(Double.parseDouble(matcher.group()));
                row_count++; // Avoid reading the last column! This will be treated next
            }
            // Process the last column separately
            String[] entryParts = entry.split("\\s+");
            String lastColumnValue = entryParts[entryParts.length - 1];

            // Add 1 for 'D' and 0 otherwise
            afields.add(lastColumnValue.contains("D") ? 1.0 : 0.0);

            // Add 1 for 'A' and 0 otherwise
            afields.add(lastColumnValue.contains("A") ? 1.0 : 0.0);
            fields.add(afields);
        }
        return fields;
    }

    public static String getCanonicalSmiles(IAtomContainer molecule) {
        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.UseAromaticSymbols);
        try {
            return smilesGenerator.create(molecule);
        } catch (CDKException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static String getMOLString(IAtomContainer molecule) throws CDKException, IOException {
        StringWriter writer = new StringWriter();
        MDLV2000Writer molWriter = new MDLV2000Writer(writer);
        molWriter.write(molecule);
        molWriter.close();
        return writer.toString();
    }

    public static String getCTABString(IAtomContainer molecule) throws CDKException, IOException {
        StringWriter writer = new StringWriter();
        SDFWriter ctabWriter = new SDFWriter(writer);
        ctabWriter.write(molecule);
        ctabWriter.close();
        return writer.toString();
    }

    public static IteratingSDFReader read(String path) {
        IteratingSDFReader reader = null;
        try {
            FileInputStream in = new FileInputStream(path);
            reader = new IteratingSDFReader(in, DefaultChemObjectBuilder.getInstance());
        } catch (Exception e) {
            System.err.println("sdf path not found: " + path);
            System.exit(1);
        }
        return reader;
    }
}