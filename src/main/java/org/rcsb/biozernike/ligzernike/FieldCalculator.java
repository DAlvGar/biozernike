package org.rcsb.biozernike.ligzernike;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.rcsb.biozernike.molecules.SDFReader;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

public class FieldCalculator {

    private static double spacing = 0.75;
    /**
     * frame
     * sets the size of the volume calculations, set a bigger frame to get more
     * space to show fields
     */
    private static final double frame = 5;
    private static double[][] exponentialFunctionParams = {
            { -0.25918178, 0.99676629 },
            { -0.19200658, 0.93042921 },
            { -0.14224198, 0.83152089 },
            { -0.10537545, 0.72138127 },
            { -0.07806405, 0.61247644 },
            { -0.05783127, 0.51156498 },
            { -0.04284246, 0.42181912 },
            { -0.03173847, 0.34422976 },
            { -0.02351244, 0.27852412 },
            { -0.01741844, 0.22375419 },
            { -0.0129039, 0.17866508 },
            { -0.00955944, 0.14191779 },
            { -0.00708181, 0.1122171 },
            { -0.00524633, 0.08837881 },
            { -0.00388658, 0.06935921 },
            { -0.00287925, 0.05426182 },
            { -0.002133, 0.04233114 },
            { -0.00158017, 0.03293985 },
            { -0.00117062, 0.02557305 },
            { -0.00086721, 0.0198122 },
            { -0.00064245, 0.01531968 },
            { -0.00047594, 0.01182504 },
            { -0.00035258, 0.00911279 },
            { -0.0002612, 0.00701212 },
            { -0.0001935, 0.00538821 },
            { -0.00014335, 0.00413503 },
            { -0.0001062, 0.0031695 },
            { -0.00007867, 0.0024267 },
            { -0.00005828, 0.00185602 },
            { -0.00004318, 0.00141815 },
            { -0.00003199, 0.00108258 },
            { -0.0000237, 0.00082569 },
            { -0.00001755, 0.00062924 },
            { -0.000013, 0.00047916 },
            { -0.00000963, 0.0003646 }
    };
    // private static CompletableFuture<String> getFieldAsync(IAtomContainer
    // molecule, int fieldID) {
    // return CompletableFuture.supplyAsync(() -> getField(molecule, fieldID));
    // }

    public static List<Volume> projectMultiField(IAtomContainer molecule, int[] multifieldID) {
        IAtomContainerExtremes moleculeExtremes = getMoleculeExtremes(molecule);
        IGridSize gridSize = getGridSize(moleculeExtremes);
        double[] corner = { gridSize.zeroX * spacing, gridSize.zeroY * spacing, gridSize.zeroZ * spacing };
        int[] dimensions = { gridSize.ix, gridSize.iy, gridSize.iz };
        List<List<Double>> atomParameters = SDFReader.getAtomFields(molecule);
        List<Volume> volumeList = new ArrayList<>();
        List<double[]> projections = getMultipleProjections(molecule, gridSize, atomParameters, multifieldID);

        for (double[] projection : projections) {
            double[] voxelArray = OpenDXIO.rowToColumnFlatten(projection, dimensions); // Transform order
            Volume volume = new Volume();
            volume.reset();
            volume.setCorner(corner);
            volume.setDimensions(dimensions);
            volume.setGridWidth(spacing);
            volume.setVoxelArray(voxelArray);
            volumeList.add(volume);
        }
        return volumeList;
    }

    public static Volume projectField(IAtomContainer molecule, int fieldID) {
        IAtomContainerExtremes moleculeExtremes = getMoleculeExtremes(molecule);
        IGridSize gridSize = getGridSize(moleculeExtremes);
        double[] corner = { gridSize.zeroX * spacing, gridSize.zeroY * spacing, gridSize.zeroZ * spacing };
        int[] dimensions = { gridSize.ix, gridSize.iy, gridSize.iz };
        List<List<Double>> atomParameters = SDFReader.getAtomFields(molecule);

        double[] projections = getProjections(molecule, gridSize, atomParameters, fieldID);
        double[] voxelArray = OpenDXIO.rowToColumnFlatten(projections, dimensions); // Transform order
        Volume volume = new Volume();
        volume.reset();
        volume.setCorner(corner);
        volume.setDimensions(dimensions);
        volume.setGridWidth(spacing);
        volume.setVoxelArray(voxelArray);
        return volume;
    }

    private static double[] getProjections(IAtomContainer molecule, IGridSize gridSize,
            List<List<Double>> atomParameters,
            int fieldID) {
        double[] projections = new double[gridSize.iNumGridPoints];
        int gridPoint = 0;

        for (int index = 0; index < gridSize.ix; index++) {
            double xcoordinate = (index + gridSize.zeroX) * spacing;

            for (int jIndex = 0; jIndex < gridSize.iy; jIndex++) {
                double ycoordinate = (jIndex + gridSize.zeroY) * spacing;

                for (int kIndex = 0; kIndex < gridSize.iz; kIndex++) {
                    double zcoordinate = (kIndex + gridSize.zeroZ) * spacing;

                    Atom atom = new Atom(xcoordinate, ycoordinate, zcoordinate);

                    projections[gridPoint] = newProjection(atom, molecule, fieldID, atomParameters) * 10; // TODO: WHY
                                                                                                          // SCALE IS
                                                                                                          // DIFFERENT??
                    gridPoint++;
                }
            }
        }

        return projections;
    }

    private static List<double[]> getMultipleProjections(IAtomContainer molecule, IGridSize gridSize,
            List<List<Double>> atomParameters,
            int[] multifieldID) {

        int N = multifieldID.length;
        List<double[]> projectionList = new ArrayList<>(N);
        for (int f : multifieldID) {
            projectionList.add(new double[gridSize.iNumGridPoints]);
        }
        int gridPoint = 0;
        double[] partial = new double[N];
        for (int index = 0; index < gridSize.ix; index++) {
            double xcoordinate = (index + gridSize.zeroX) * spacing;

            for (int jIndex = 0; jIndex < gridSize.iy; jIndex++) {
                double ycoordinate = (jIndex + gridSize.zeroY) * spacing;

                for (int kIndex = 0; kIndex < gridSize.iz; kIndex++) {
                    double zcoordinate = (kIndex + gridSize.zeroZ) * spacing;

                    Atom atom = new Atom(xcoordinate, ycoordinate, zcoordinate);

                    partial = newProjectionMulti(atom, molecule, multifieldID, atomParameters); // TODO:

                    for (int k = 0; k < N; k++) {
                        projectionList.get(k)[gridPoint] = partial[k] * 10;
                    }

                    gridPoint++;
                }
            }
        }

        return projectionList;
    }

    private static IGridSize getGridSize(IAtomContainerExtremes moleculeExtremes) {

        IGridSize gridSize = new IGridSize();
        gridSize.xmax = moleculeExtremes.xmax + frame;
        gridSize.xmin = moleculeExtremes.xmin - frame;
        gridSize.ymax = moleculeExtremes.ymax + frame;
        gridSize.ymin = moleculeExtremes.ymin - frame;
        gridSize.zmax = moleculeExtremes.zmax + frame;
        gridSize.zmin = moleculeExtremes.zmin - frame;

        gridSize.zeroX = (int) Math.ceil(gridSize.xmin / spacing);
        gridSize.zeroY = (int) Math.ceil(gridSize.ymin / spacing);
        gridSize.zeroZ = (int) Math.ceil(gridSize.zmin / spacing);

        gridSize.ix = (int) Math.ceil((gridSize.xmax / spacing) - gridSize.zeroX + 1);
        gridSize.iy = (int) Math.ceil((gridSize.ymax / spacing) - gridSize.zeroY + 1);
        gridSize.iz = (int) Math.ceil((gridSize.zmax / spacing) - gridSize.zeroZ + 1);

        gridSize.iNumGridPoints = gridSize.ix * gridSize.iy * gridSize.iz;

        return gridSize;
    }

    private static IAtomContainerExtremes getMoleculeExtremes(IAtomContainer molecule) {
        IAtomContainerExtremes extremes = new IAtomContainerExtremes();
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            double x = molecule.getAtom(i).getPoint3d().x;
            double y = molecule.getAtom(i).getPoint3d().y;
            double z = molecule.getAtom(i).getPoint3d().z;

            if (x > extremes.xmax) {
                extremes.xmax = x;
            }
            if (x < extremes.xmin) {
                extremes.xmin = x;
            }
            if (y > extremes.ymax) {
                extremes.ymax = y;
            }
            if (y < extremes.ymin) {
                extremes.ymin = y;
            }
            if (z > extremes.zmax) {
                extremes.zmax = z;
            }
            if (z < extremes.zmin) {
                extremes.zmin = z;
            }
        }
        return extremes;
    }

    private static double[] newProjectionMulti(Atom atom, IAtomContainer molecule, int[] multifieldID,
            List<List<Double>> atomParameters) {
        double sqrDistance;
        double[] projection = new double[multifieldID.length];
        double minimumDist = 999;
        Atom newAtom;

        for (int lIndex = 0; lIndex < molecule.getAtomCount(); lIndex++) {
            newAtom = new Atom(molecule.getAtom(lIndex).getPoint3d().x,
                    molecule.getAtom(lIndex).getPoint3d().y,
                    molecule.getAtom(lIndex).getPoint3d().z);
            sqrDistance = atom.squareEuclideanDistance(newAtom);

            for (int i = 0; i < multifieldID.length; i++) {
                int fieldID = multifieldID[i];
                double param = atomParameters.get(lIndex).get(fieldID);
                if (fieldID < 5) { // below 5 is charge and hydrophobic numbers, columns 6th and 7th are hbonds
                    projection[i] += getExponentialProjection(param, sqrDistance);
                } else if (sqrDistance < minimumDist
                        && sqrDistance < 1.5
                        && Math.abs(param) == 1) {
                    minimumDist = sqrDistance;
                    projection[i] = param / sqrDistance;
                }
            }
        }

        return projection;
    }

    private static double newProjection(Atom atom, IAtomContainer molecule, int fieldID,
            List<List<Double>> atomParameters) {
        double sqrDistance;
        double projection = 0.0;
        double minimumDist = 999;
        Atom newAtom;

        for (int lIndex = 0; lIndex < molecule.getAtomCount(); lIndex++) {
            newAtom = new Atom(molecule.getAtom(lIndex).getPoint3d().x,
                    molecule.getAtom(lIndex).getPoint3d().y,
                    molecule.getAtom(lIndex).getPoint3d().z);
            sqrDistance = atom.squareEuclideanDistance(newAtom);

            double param = atomParameters.get(lIndex).get(fieldID);
            if (fieldID < 5) { // below 5 is charge and hydrophobic numbers, columns 6th and 7th are hbonds
                projection += getExponentialProjection(param, sqrDistance);
            } else if (sqrDistance < minimumDist
                    && sqrDistance < 1.5
                    && Math.abs(param) == 1) {
                minimumDist = sqrDistance;
                projection = param / sqrDistance;
            }
        }

        return projection;
    }

    private static String fieldFileFormat(IGridSize gridSize, double[] projections) {
        StringBuilder DXFile = new StringBuilder();
        double originX = gridSize.zeroX * spacing;
        double originY = gridSize.zeroY * spacing;
        double originZ = gridSize.zeroZ * spacing;

        DXFile.append("# Data from 1.4.1\n");
        DXFile.append("# \n");
        DXFile.append("# POTENTIAL (kT/e)\n");
        DXFile.append("# \n");

        DXFile.append("object 1 class gridpositions counts ")
                .append(gridSize.ix).append(" ")
                .append(gridSize.iy).append(" ")
                .append(gridSize.iz).append("\n");

        DXFile.append("origin ").append(originX).append(" ").append(originY).append(" ").append(originZ).append("\n");
        DXFile.append("delta ").append(spacing).append(" 0.000000e+00 0.000000e+00\n");
        DXFile.append("delta 0.000000e+00 ").append(spacing).append(" 0.000000e+00\n");
        DXFile.append("delta 0.000000e+00 0.000000e+00 ").append(spacing).append("\n");
        DXFile.append("object 2 class gridconnections counts ")
                .append(gridSize.ix).append(" ")
                .append(gridSize.iy).append(" ")
                .append(gridSize.iz).append("\n");

        DXFile.append("object 3 class array type double rank 0 items ")
                .append(gridSize.iNumGridPoints).append(" data follows\n");

        int gridPoint = 0;
        for (int index = 0; index < gridSize.iNumGridPoints / 3; index++) {
            for (int jindex = 0; jindex < 3; jindex++) {
                DXFile.append(projections[gridPoint]).append(" ");
                gridPoint++;
            }
            DXFile.append("\n");
        }

        for (int index = 0; index < gridSize.iNumGridPoints % 3; index++) {
            DXFile.append(projections[gridPoint]).append(" ");
            gridPoint++;
        }

        return DXFile.toString();
    }

    private static double getExponentialProjection(double atomParameter, double sqrDistance) {
        int exponentialIndex = (int) Math.round(sqrDistance);
        if (exponentialIndex < 0) {
            throw new RuntimeException("Exponential function error");
        }

        if (exponentialIndex >= exponentialFunctionParams.length) {
            // When more than the available parameters, assume the projection is 0
            return 0.0;
        }

        return (exponentialFunctionParams[exponentialIndex][1] * atomParameter)
                + (exponentialFunctionParams[exponentialIndex][0] * sqrDistance * atomParameter);
    }

    private static class Atom {
        private final double x;
        private final double y;
        private final double z;

        public Atom(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public double squareEuclideanDistance(Atom other) {
            return Math.pow(x - other.x, 2) + Math.pow(y - other.y, 2) + Math.pow(z - other.z, 2);
        }
    }

    private static class IGridSize {
        private double xmax;
        private double xmin;
        private double ymax;
        private double ymin;
        private double zmax;
        private double zmin;
        private int zeroX;
        private int zeroY;
        private int zeroZ;
        private int ix;
        private int iy;
        private int iz;
        private int iNumGridPoints;
    }

    private static class IAtomContainerExtremes {
        private double xmax;
        private double xmin;
        private double ymax;
        private double ymin;
        private double zmax;
        private double zmin;
    }
}
