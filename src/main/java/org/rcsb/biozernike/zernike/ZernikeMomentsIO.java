package org.rcsb.biozernike.zernike;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.rcsb.biozernike.complex.Complex;

public class ZernikeMomentsIO {

    public static void writeDouble(String filePath, List<Double> doubleNumbers) throws IOException{
        BufferedWriter writer = new BufferedWriter(new FileWriter(filePath));
        for (Double number : doubleNumbers) {
            writer.write(number.toString() + System.lineSeparator());
        }
        writer.close();
    }

    public static List<Double> readDouble(String filePath) throws IOException, ClassNotFoundException{
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

    public static List<Complex> readComplex(String filePath, boolean binary) throws ClassNotFoundException{
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
                    // Assuming the format is "(real + imaginary i)"
                    String[] parts = line.split("[+i()]");
                    double real = Double.parseDouble(parts[1].trim());
                    double imaginary = Double.parseDouble(parts[2].trim());
                    complexNumbers.add(new Complex(real, imaginary));
                }
                reader.close();
            }
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
        return complexNumbers;
    }

}
