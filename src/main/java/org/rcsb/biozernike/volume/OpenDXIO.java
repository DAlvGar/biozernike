package org.rcsb.biozernike.volume;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

/**
 * OpenDXParser
 */
public class OpenDXIO {
	/**
	 * Read volume from given file in specified format
	 * 
	 * @param file the file
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	public static Volume read(String fileName) throws IOException {
		return read_dx(new File(fileName));
	}

	/**
	 * Read volume from given input stream in DX format
	 * 
	 * @param file the file
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	private static Volume read_dx(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;

		double[] origin = null;
		int[] dims = new int[3];
		double[] voxels = null;
		double gridWidth = 0.0;

		int valuesIndex = 0;
		boolean readingValues = false;

		while ((line = reader.readLine()) != null) {
			if (!readingValues) {
				if (line.startsWith("object") && line.contains("gridpositions")) {
					// Parse gridpositions
					String[] parts = line.split(" ");
					for (int i = 0; i < 3; i++) {
						dims[i] = Integer.parseInt(parts[i + 5]);
					}
				} else if (line.contains("origin")) {
					// Parse origin
					String[] parts = line.split(" ");
					origin = new double[parts.length - 1];
					for (int i = 0; i < origin.length; i++) {
						origin[i] = Float.parseFloat(parts[i + 1]);
					}
				} else if (line.contains("delta")) {
					// Parse spacing (delta)
					String[] parts = line.split(" +");
					for (int i = 0; i < parts.length - 1; i++) {
						double width = Double.parseDouble(parts[i + 1]);
						if (width > gridWidth) {
							gridWidth = width;
						}
					}
				} else if (line.contains("array")) {
					// Initialize values array based on dimensions
					String[] parts = line.split(" +");
					int totalValues = Integer.parseInt(parts[9]);
					voxels = new double[totalValues];
					readingValues = true;
				}
			} else {
				String[] parts = line.trim().split("\\s+");
				for (String part : parts) {
					voxels[valuesIndex++] = Double.parseDouble(part);
				}
			}
		}
		reader.close();
		// Transform row-first dx format to column-first ccp4 format used by the rest of the program
		double[] voxels_col = rowToColumnFlatten(voxels, dims);
		Volume volume = new Volume();
		volume.createFromData(dims, voxels_col, gridWidth);
		volume.setCorner(origin);
		return volume;
	}

	/**
	 * Write DX file from given input Volume
	 * 
	 * @param file          output file
	 * @param volume        Volume to write
	 * @param origin        grid origin position
	 * @param delta         grid spacing in each direction
	 * @param valuesPerLine Usually 3
	 * @return void
	 * @throws IOException if problems reading file
	 */
	public static void write(String fileName, Volume volume) throws IOException {
		File file = new File(fileName);
		int valuesPerLine = 3;
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		double delta = volume.getGridWidth();
		double[] delta1 = { delta, 0, 0 };
		double[] delta2 = { 0, delta, 0 };
		double[] delta3 = { 0, 0, delta };
		// Writing header information including dimensions, origin, and delta
		// Transform colum-first voxel array to row-first in DX format
		double[] voxels = columnToRowFlatten(volume.getVoxelArray(),volume.getDimensions());
		writer.write("object 1 class gridpositions counts " + volume.getDimensions()[0] + " " +
				volume.getDimensions()[1] + " " + volume.getDimensions()[2] + "\n");
		writer.write("origin " + formatArray(volume.getCorner()) + "\n");
		writer.write("delta " + formatArray(delta1) + "\n");
		writer.write("delta " + formatArray(delta2) + "\n");
		writer.write("delta " + formatArray(delta3) + "\n");
		writer.write("object 2 class gridconnections counts " + volume.getDimensions()[0] + " " +
				volume.getDimensions()[1] + " " + volume.getDimensions()[2] + "\n");
		writer.write("object 3 class array type double rank 0 items " + voxels.length + " data follows\n");

		// Writing the voxel values
		for (int i = 0; i < voxels.length; i += valuesPerLine) {
			int endIndex = Math.min(i + valuesPerLine, voxels.length);
			for (int j = i; j < endIndex; j++) {
				writer.write(String.format(Locale.US, "%.3f", voxels[j]) + " ");
			}
			writer.write("\n");
		}

		writer.close();
	}

	// Helper method to format arrays for output
	private static String formatArray(double[] array) {
		StringBuilder builder = new StringBuilder();
		for (double value : array) {
			builder.append(String.format(Locale.US, "%.2f", value)).append(" ");
		}
		return builder.toString().trim();
	}

	public static double[] rowToColumnFlatten(double[] rowFlatten, int[] dims) {
		if (dims.length != 3) {
			throw new IllegalArgumentException("Expected a 3D shape");
		}

		if (rowFlatten.length != dims[0] * dims[1] * dims[2]) {
			throw new IllegalArgumentException("Input array size doesn't match the given shape");
		}

		int rows = dims[0];
		int cols = dims[1];
		int depth = dims[2];

		double[] columnFlatten = new double[rowFlatten.length];

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < depth; k++) {
					int indexColumn = k * rows * cols + j * rows + i;
					int indexRow = i * cols * depth + j * depth + k;
					columnFlatten[indexColumn] = rowFlatten[indexRow];
				}
			}
		}
		return columnFlatten;
	}


	private static double[] columnToRowFlatten(double[] columnFlatten, int[] dims) {
		if (dims.length != 3) {
			throw new IllegalArgumentException("Expected a 3D shape");
		}

		if (columnFlatten.length != dims[0] * dims[1] * dims[2]) {
			throw new IllegalArgumentException("Input array size doesn't match the given shape");
		}

		int rows = dims[0];
		int cols = dims[1];
		int depth = dims[2];

		double[] rowFlatten = new double[columnFlatten.length];

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < depth; k++) {
					int indexColumn = k * rows * cols + j * rows + i;
					int indexRow = i * cols * depth + j * depth + k;
					rowFlatten[indexRow] = columnFlatten[indexColumn];
				}
			}
		}
		return rowFlatten;
	}


}