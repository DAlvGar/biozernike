package org.rcsb.biozernike.zernike;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.rcsb.biozernike.complex.Complex;

import java.util.regex.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class ZernikeMomentsIO {

    private static String regexPattern = "(-?\\d+(\\.\\d+)?(?:E[+-]?\\d+)?)([+-]\\d+(\\.\\d+)?(?:E[+-]?\\d+)?)j";
    private static Pattern pattern = Pattern.compile(regexPattern);

    private static Complex strToComplex(String s) {
        Matcher matcher = pattern.matcher(s);
        if (matcher.find()) {
            double real = Double.parseDouble(matcher.group(1));
            double imaginary = Double.parseDouble(matcher.group(3));
            return new Complex(real, imaginary);
        } else {
            throw new IllegalArgumentException("Invalid complex number format: " + s);
        }
    }

    public static void writeDouble(String filePath, List<Double> doubleNumbers) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(filePath));
        for (Double number : doubleNumbers) {
            writer.write(number.toString() + System.lineSeparator());
        }
        writer.close();
    }

    public static String complexListToString(List<Complex> complex){
        return complex.stream().map(String::valueOf).collect(Collectors.joining(","));
    }

    public static List<Double> readDouble(String filePath) throws IOException, ClassNotFoundException {
        List<Double> doubleNumbers = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                // Assuming the format is "(real + imaginary i)"
                double d = Double.parseDouble(line.trim());
                doubleNumbers.add(d);
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return doubleNumbers;
    }

    public static List<Double> readDoubleFromString(String doubleStr){
        List<Double> doubleList = Arrays.asList(Stream.of(doubleStr.split(",")).mapToDouble(Double::parseDouble).boxed().toArray(Double[]::new));
        return doubleList;
    }
    
    public static void writeComplex(String filePath, List<Complex> complexNumbers, boolean binary) throws IOException {
        try {
            if (binary) {
                ObjectOutputStream outputStream = new ObjectOutputStream(new FileOutputStream(filePath));
                outputStream.writeObject(complexNumbers);
                outputStream.close();
            } else {
                BufferedWriter writer = new BufferedWriter(new FileWriter(filePath));
                for (Complex number : complexNumbers) {
                    writer.write(number.toString() + System.lineSeparator());
                }
                writer.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<Complex> readComplex(String filePath, boolean binary) throws ClassNotFoundException {
        List<Complex> complexNumbers = new ArrayList<>();
        try {
            if (binary) {
                ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(filePath));
                Object obj = inputStream.readObject();
                if (obj instanceof List) {
                    complexNumbers = (List<Complex>) obj;
                }
                inputStream.close();
            } else {
                BufferedReader reader = new BufferedReader(new FileReader(filePath));
                String line;
                while ((line = reader.readLine()) != null) {
                    complexNumbers.add(strToComplex(line));
                }
                reader.close();
            }
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
        return complexNumbers;
    }

    public static List<Complex> readComplex(String concatenatedComplex) {
        List<Complex> complexNumbers = new ArrayList<>();

        String[] splitComplex = concatenatedComplex.split(",");
        String line;
        for (int i = 0; i < splitComplex.length; i++) {
            line = splitComplex[i];
            complexNumbers.add(strToComplex(line));
        }
        return complexNumbers;
    }

}
