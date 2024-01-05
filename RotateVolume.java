package org.rcsb.biozernike.volume;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.TrivariateGridInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.stream.IntStream;

import javax.vecmath.Matrix4d;

public class RotateVolume {

    public static void rotate(String inputFilename, String outputFilename, Matrix4d transformationMatrix) {

        try {
            Volume originalVolume = OpenDXIO.read(inputFilename);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public static void rotate(Volume originalVolume, String outputFilename, Matrix4d transformationMatrix) {

        int[] dimersions = originalVolume.getDimensions();
        // Create spatial coordinates based on the grid dimensions
        double[] x = IntStream.range(0, dimersions[0]).asDoubleStream().toArray();
        double[] y = IntStream.range(0, dimersions[1]).asDoubleStream().toArray();
        double[] z = IntStream.range(0, dimersions[2]).asDoubleStream().toArray();

        // Rotate the coordinates
        Rotation rotation = new Rotation(transformationMatrix, RotationConvention.VECTOR_OPERATOR);
        double[][] rotatedCoordinates = new double[x.length * y.length * z.length][3];
        int idx = 0;
        for (double xi : x) {
            for (double yi : y) {
                for (double zi : z) {
                    Vector3D v = rotation.applyTo(new Vector3D(xi, yi, zi));
                    rotatedCoordinates[idx++] = v.toArray();
                }
            }
        }

            // Interpolate the original density map onto the new rotated coordinates
            TrivariateGridInterpolator interpolator = new LinearInterpolator();
            PolynomialSplineFunction[][][] interpolatedValues = interpolator.interpolate(new double[][][]{x, y, z}, originalDensityMap);

            double[][][] newDensityMap = new double[gridDimensions[0]][gridDimensions[1]][gridDimensions[2]];
            for (int i = 0; i < rotatedCoordinates.length; i++) {
                double[] coord = rotatedCoordinates[i];
                double val = interpolatedValues[0][0][0].value(coord);
                int xIdx = i / (gridDimensions[1] * gridDimensions[2]);
                int yIdx = (i / gridDimensions[2]) % gridDimensions[1];
                int zIdx = i % gridDimensions[2];
                newDensityMap[xIdx][yIdx][zIdx] = val;
            }

            // Define origin and delta
            double[] origin = {0, 0, 0};
            double delta = 1.0;

            // Write the rotated density map to a new DX file
            writeDXFile(outputFilename, newDensityMap, gridDimensions, origin, delta);
        
    }
}
