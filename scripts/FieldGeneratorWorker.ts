import Atom from '../model/Atom';
import type IGridSize from '@/model/IGridSize';
import type { IMoleculeExtremes } from '@/model/IMoleculeExtremes';
import * as OCL from 'openchemlib';
import Constants from '../plugins/Constants';

/**
 * spacing
 * @type {number}
 * sets the distance between point to calculate the projections
 * ! Be careful, with smaller spacing the computation time grows exponentially
 */
let spacing: number = 0.75; // ? Factor de correcció? Amstrongs?
// let spacing: number = 0.15; // ? Factor de correcció? Amstrongs?
/**
 * frame
 * @type {number}
 * sets the size of the volume calculations, set a bigger frame to get more space to show fields
 * ? Marks the grid size anything outside of it never will be rendered
 */
const frame: number = 5;

// Linearization of exponential function with R^2 = 0.9981
const exponentialFunctionParams: Array<Array<number>>
= [
	[ -0.25918178, 0.99676629 ],
	[ -0.19200658, 0.93042921 ],
	[ -0.14224198, 0.83152089 ],
	[ -0.10537545, 0.72138127 ],
	[ -0.07806405, 0.61247644 ],
	[ -0.05783127, 0.51156498 ],
	[ -0.04284246, 0.42181912 ],
	[ -0.03173847, 0.34422976 ],
	[ -0.02351244, 0.27852412 ],
	[ -0.01741844, 0.22375419 ],
	[ -0.0129039, 0.17866508 ],
	[ -0.00955944, 0.14191779 ],
	[ -0.00708181, 0.1122171 ],
	[ -0.00524633, 0.08837881 ],
	[ -0.00388658, 0.06935921 ],
	[ -0.00287925, 0.05426182 ],
	[ -0.002133, 0.04233114 ],
	[ -0.00158017, 0.03293985 ],
	[ -0.00117062, 0.02557305 ],
	[ -0.00086721, 0.0198122 ],
	[ -0.00064245, 0.01531968 ],
	[ -0.00047594, 0.01182504 ],
	[ -0.00035258, 0.00911279 ],
	[ -0.0002612, 0.00701212 ],
	[ -0.0001935, 0.00538821 ],
	[ -0.00014335, 0.00413503 ],
	[ -0.0001062, 0.0031695 ],
	[ -0.00007867, 0.0024267 ],
	[ -0.00005828, 0.00185602 ],
	[ -0.00004318, 0.00141815 ],
	[ -0.00003199, 0.00108258 ],
	[ -0.0000237, 0.00082569 ],
	[ -0.00001755, 0.00062924 ],
	[ -0.000013, 0.00047916 ],
	[ -0.00000963, 0.0003646 ],
];

self.onmessage = function (event) {
	const dataVal = JSON.parse(event.data);

	_getField(
		dataVal.moleculeFile,
		dataVal.atomParameters,
		dataVal.iFieldID,
	)
		.then(dxFile => {
			if (dxFile !== '') {
				self.postMessage(dxFile);
			}
		});
};

/**
 * getField Parse the moleculeFile into a sdx field file in order to show a field
 * @param moleculeFile {string} A string formated with the sdf molecule file structure
 * @param atomParameters {string[][]} Atoms parameters parsed from the molecule file
 * @param fieldID {number} The file type selected
 * @returns {Promise<string>} Returns a promisa of a string
 */
const _getField = function (
	moleculeFile: string,
	atomParameters: Array<Array<string>>,
	fieldID: number,
): Promise<string> {
	// registerPromiseWorker.
	const promise: Promise<string> = new Promise((resolve, reject) => {
		if (!moleculeFile) {
			const err = new Error('No molecule file selected');

			reject(err);
			throw err;
		}

		// Parse the sdf molecule file type to an sd field type file
		const parser = new OCL.SDFileParser(moleculeFile, []);

		if (fieldID < 0 || fieldID >= Constants.MAX_NUMBER_OF_FIELDS) {
			const err = new Error('Unknown field type: ' + fieldID);

			reject(err);
			throw err;
		}

		if (atomParameters && atomParameters.length > 0 && fieldID >= atomParameters[0].length) {
			const err = new Error('Data not available for selected field (' + fieldID + ')');

			reject(err);
			throw err;
		}

		if (fieldID >= Constants.NUMBER_OF_FIELDS) {
			spacing = 0.25;
		}

		let molecule: OCL.Molecule;
		let moleculeExtemes: IMoleculeExtremes;
		let gridSize: IGridSize = {
			xmax: 0,
			xmin: 0,
			ymax: 0,
			ymin: 0,
			zmax: 0,
			zmin: 0,
			zeroX: 0,
			zeroY: 0,
			zeroZ: 0,
			ix: 0,
			iy: 0,
			iz: 0,
			iNumGridPoints: 0,
		};

		let projections: Array<number> = [];

		// Get the projections of the molecule
		while (parser.next()) {
			molecule = parser.getMolecule();
			moleculeExtemes = _getMoleculeExtremes(molecule);
			gridSize = _getGridSize(moleculeExtemes);
			projections = _getProjections(
				molecule,
				gridSize,
				atomParameters,
				fieldID,
			);
		}

		const dXFile: string = _fieldFileFormat(
			gridSize,
			projections,
		);

		resolve(dXFile);
	});

	return promise;
};

/**
 * Molecules get molecule furthests points, in oder to know the limits of the molecule
 * @param molecule {OCL.Molecule} The molecule that is to be mesured
 * @returns
 */
const _getMoleculeExtremes = function (molecule: OCL.Molecule): IMoleculeExtremes {
	const extremes: IMoleculeExtremes
	= {
		xmax: molecule.getAtomX(0),
		xmin: molecule.getAtomX(0),
		ymax: -molecule.getAtomY(0),
		ymin: -molecule.getAtomY(0),
		zmax: -molecule.getAtomZ(0),
		zmin: -molecule.getAtomZ(0),
	};

	for (let index = 0; index < molecule.getAllAtoms(); index++) {
		if (molecule.getAtomX(index) > extremes.xmax) {
			extremes.xmax = molecule.getAtomX(index);
		}

		if (molecule.getAtomX(index) < extremes.xmin) {
			extremes.xmin = molecule.getAtomX(index);
		}

		if (-molecule.getAtomY(index) > extremes.ymax) {
			extremes.ymax = -molecule.getAtomY(index);
		}

		if (-molecule.getAtomY(index) < extremes.ymin) {
			extremes.ymin = -molecule.getAtomY(index);
		}

		if (-molecule.getAtomZ(index) > extremes.zmax) {
			extremes.zmax = -molecule.getAtomZ(index);
		}

		if (-molecule.getAtomZ(index) < extremes.zmin) {
			extremes.zmin = -molecule.getAtomZ(index);
		}
	}

	//! WARNING: Atom coordinates from Y and Z seem to be with opposite sign in this lib
	return extremes;
};

/**
 * Extremes get grid size, Calculate the frame size, using the size of the actual molecule
 * @param extremes {IMolecueExtremes} The molecule size
 * @returns
 */
const _getGridSize = function (extremes: IMoleculeExtremes): IGridSize {
	// Calculating the frame arround which the field will be calculated
	const gridSize: IGridSize
	= {
		xmax: extremes.xmax + frame,
		xmin: extremes.xmin - frame,
		ymax: extremes.ymax + frame,
		ymin: extremes.ymin - frame,
		zmax: extremes.zmax + frame,
		zmin: extremes.zmin - frame,
		zeroX: 0,
		zeroY: 0,
		zeroZ: 0,
		ix: 0,
		iy: 0,
		iz: 0,
		iNumGridPoints: 0,
	};

	// Zero is calculated to ensure grid points of different molecules are aligned
	gridSize.zeroX = Math.ceil(gridSize.xmin / spacing);
	gridSize.zeroY = Math.ceil(gridSize.ymin / spacing);
	gridSize.zeroZ = Math.ceil(gridSize.zmin / spacing);

	// Calculating number of elements in the grid in every coordinate
	// max / spacing - min / spacing + 1
	gridSize.ix = Math.ceil((gridSize.xmax / spacing) - gridSize.zeroX + 1);
	gridSize.iy = Math.ceil((gridSize.ymax / spacing) - gridSize.zeroY + 1);
	gridSize.iz = Math.ceil((gridSize.zmax / spacing) - gridSize.zeroZ + 1);

	// Calculating actual size of the grid
	gridSize.iNumGridPoints = gridSize.ix * gridSize.iy * gridSize.iz;

	return gridSize;
};

/**
 * Molecules get projections, calculates all the projections of the atoms in a mesh
 * @param molecule {OCL.Molecule}
 * @param gridSize {IGridSize}
 * @param atomParameters {string[][]}
 * @param fieldID {number}
 * @returns
 */
const _getProjections = function (
	molecule: OCL.Molecule,
	gridSize: IGridSize,
	atomParameters: Array<Array<string>>,
	fieldID: number,
): Array<number> {
	let xcoordinate: number;
	let ycoordinate: number;
	let zcoordinate: number;
	const projections: Array<number> = [];
	let atom: Atom;

	// Creating the main loop of the gird
	for (let index = 0; index < gridSize.ix; index++) {
		xcoordinate = (index + gridSize.zeroX) * spacing; // The position multiplied

		for (let jIndex = 0; jIndex < gridSize.iy; jIndex++) {
			ycoordinate = (jIndex + gridSize.zeroY) * spacing;

			for (let kIndex = 0; kIndex < gridSize.iz; kIndex++) {
				zcoordinate = (kIndex + gridSize.zeroZ) * spacing;

				atom = new Atom(xcoordinate, ycoordinate, zcoordinate);

				projections.push(
					_newProjection(
						molecule,
						atom,
						fieldID,
						atomParameters,
					),
				);
			}
		}
	}

	return projections;
};

/**
 * Molecules new projection, calculates the projection using the exponential projection
 * @param molecule {OCL.Molecule}
 * @param atom {Atom}
 * @param fieldID {number}
 * @param atomParameters {string[][]}
 * @returns
 */
const _newProjection = function (
	molecule: OCL.Molecule,
	atom: Atom,
	fieldID: number,
	atomParameters: Array<Array<string>>,
): number {
	let sqrDistance: number;
	let projection: number = 0.0;
	let minimumDist: number = 999;
	let newAtom: Atom;

	for (let lIndex = 0; lIndex < molecule.getAllAtoms(); lIndex++) {
		newAtom = new Atom(
			molecule.getAtomX(lIndex),
			-molecule.getAtomY(lIndex),
			-molecule.getAtomZ(lIndex),
		);
		sqrDistance = atom.squareEuclideanDistance(newAtom);
		if (fieldID < Constants.NUMBER_OF_FIELDS) {
			projection += _getExponentialProjection(
				parseFloat(atomParameters[lIndex][fieldID]),
				sqrDistance,
			);
		} else if (
			sqrDistance < minimumDist
			&& sqrDistance < 1.5
			&& Math.abs(parseFloat(atomParameters[lIndex][fieldID])) === 1
		) {
			minimumDist = sqrDistance;
			projection = parseFloat(atomParameters[lIndex][fieldID]) / sqrDistance;
		}
	}

	return projection;
};

const _fieldFileFormat = function (
	gridSize: IGridSize,
	projections: Array<number>,
) {
	const originX: number = gridSize.zeroX * spacing;
	const originY: number = gridSize.zeroY * spacing;
	const originZ: number = gridSize.zeroZ * spacing;

	let DXFile: string;
	DXFile = '# Data from 1.4.1\n';
	DXFile += '# \n';
	DXFile += '# POTENTIAL (kT/e)\n';
	DXFile += '# \n';

	DXFile
		+= 'object 1 class gridpositions counts '
		+ gridSize.ix
		+ ' '
		+ gridSize.iy
		+ ' '
		+ gridSize.iz
		+ '\n';

	DXFile += 'origin ' + originX + ' ' + originY + ' ' + originZ + '\n';
	DXFile += 'delta ' + spacing + ' 0.000000e+00 0.000000e+00\n';
	DXFile += 'delta 0.000000e+00 ' + spacing + ' 0.000000e+00\n';
	DXFile += 'delta 0.000000e+00 0.000000e+00 ' + spacing + '\n';
	DXFile
		+= 'object 2 class gridconnections counts '
		+ gridSize.ix
		+ ' '
		+ gridSize.iy
		+ ' '
		+ gridSize.iz
		+ '\n';

	DXFile
        += 'object 3 class array type double rank 0 items '
        + gridSize.iNumGridPoints
        + ' data follows\n';

	let gridPoint: number = 0;
	for (let index = 0; index < Math.floor(gridSize.iNumGridPoints / 3); index++) {
		// DXFile += '9.000000e-00 9.000000e-00 9.000000e-00 \n';
		for (let jindex = 0; jindex < 3; jindex++) {
			DXFile += projections[gridPoint] + ' ';
			gridPoint++;
		}

		DXFile += '\n';
	}

	for (let index = 0; index < gridSize.iNumGridPoints % 3; index++) {
		// DXFile += '9.000000e-00 ';
		DXFile += projections[gridPoint] + ' ';
		gridPoint++;
	}

	return DXFile;
};

/**
 * Set the projection using an exponential function table position
 * @param atomParameter {number} atom information
 * @param sqrDistance {number} distance between the actual atom and the new atom
 * @returns
 */
const _getExponentialProjection = function (atomParameter: number, sqrDistance: number): number {
	// let exponentialIndex = Math.floor(sqrDistance); // ! Originally is using Math.floor
	const exponentialIndex = Math.round(sqrDistance);

	if (exponentialIndex < 0) {
		throw new Error('Exponential function error');
	}

	if (exponentialIndex >= 35) {
		// When more than 25A (25 x 10 e -10) we assume the projection is 0
		return 0.0;
	}

	return ((exponentialFunctionParams[exponentialIndex][1] * atomParameter)
			+ (exponentialFunctionParams[exponentialIndex][0] * sqrDistance * atomParameter));
};
