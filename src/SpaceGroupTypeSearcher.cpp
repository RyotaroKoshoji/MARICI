#include "SpaceGroupTypeSearcher.h"

#include <iostream>
#include <algorithm>
#include <map>

#include "spglib.h"

#include "ApplicationException.h"

#include "ElementSymbol.h"

using namespace ChemToolkit::Crystallography::Symmetry;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

double SpaceGroupTypeSearcher::s_defaultPrecision{ ChemToolkit::Crystallography::Symmetry::CrystallographicSymmetryOperation::defaultPrecision() };



SpaceGroupTypeSearcher::SpaceGroupTypeSearcher() noexcept
	: _precision{ s_defaultPrecision }
{
}

SpaceGroupTypeSearcher::SpaceGroupTypeSearcher(const double val)
	: _precision{ val }
{
	if (_precision <= 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), " constructor", "Input precision is not more than zero." };
}

// Constructors
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private methods

ChemToolkit::Crystallography::UnitCell SpaceGroupTypeSearcher::getDelaunayReducedUnitCell(const UnitCell& cell) const
{
	double lattice[3][3];
	{
		lattice[0][0] = cell.basisVectors()(0, 0);
		lattice[0][1] = cell.basisVectors()(0, 1);
		lattice[0][2] = cell.basisVectors()(0, 2);
		lattice[1][0] = cell.basisVectors()(1, 0);
		lattice[1][1] = cell.basisVectors()(1, 1);
		lattice[1][2] = cell.basisVectors()(1, 2);
		lattice[2][0] = cell.basisVectors()(2, 0);
		lattice[2][1] = cell.basisVectors()(2, 1);
		lattice[2][2] = cell.basisVectors()(2, 2);
	}


	if (spg_delaunay_reduce(lattice, _precision) == 0)
		throw System::ExceptionServices::ApplicationException{ typeid(*this), "getReducedUnitCell", "Could not apply delaunay reduction." };

	else
	{
		UnitCell reducedUnitCell;
		reducedUnitCell.basisVectors() = toNumericalMatrix(lattice);

		return reducedUnitCell;
	}
}

SpaceGroupTypeSearcher::SpglibOutputVariables SpaceGroupTypeSearcher::getPrimitiveOutputVariables(const SpglibInputVariables& inputVariables) const
{
	SpglibOutputVariables outputVariables = getConventionalOutputVariables(inputVariables);
	{
		SpglibInputVariables inputVariables = toSpglibInputVariables(outputVariables);
		int numAtoms = spg_standardize_cell(inputVariables.lattice, inputVariables.positions.get(), inputVariables.types.get(), inputVariables.numAtoms, 1, 1, s_defaultPrecision);


		if (numAtoms == 0)
			throw System::ExceptionServices::ApplicationException{ typeid(*this), "getPrimitiveInputVariables", "Could not find primitive structure." };

		else
		{
			const UnitCell prevUnitCell = outputVariables.unitCell;

			outputVariables.unitCell.basisVectors() = toNumericalMatrix(inputVariables.lattice);
			const NumericalMatrix inverseBasisVectors = outputVariables.unitCell.getInverseBasisVectors();
			{
				std::vector<AtomicNumber> atomicNumbers;
				std::vector<NumericalVector> fractionalCoordinates;
				std::vector<std::string> siteLabels;
				std::vector<WyckoffSymbol> wyckoffSymbols;
				std::vector<SiteSymmetrySymbol> siteSymmetrySymbols;
				std::vector<size_type> mappingToInputIndices;
				{
					constexpr double samePositionDistanceRange = 0.001;


					for (size_type index = 0; index < outputVariables.atomicNumbers.size(); ++index)
					{
						NumericalVector originalFractionalCoordinate = inverseBasisVectors * (prevUnitCell.basisVectors() * outputVariables.fractionalCoordinates.at(index));
						{
							int floorValueA = static_cast<int>(std::floor(originalFractionalCoordinate[0]));
							int floorValueB = static_cast<int>(std::floor(originalFractionalCoordinate[1]));
							int floorValueC = static_cast<int>(std::floor(originalFractionalCoordinate[2]));

							if (floorValueA != 0)
								originalFractionalCoordinate[0] -= static_cast<double>(floorValueA);
							if (floorValueB != 0)
								originalFractionalCoordinate[1] -= static_cast<double>(floorValueB);
							if (floorValueC != 0)
								originalFractionalCoordinate[2] -= static_cast<double>(floorValueC);
						}


						bool hasSameFractionalCoordinate = false;
						{
							for (int indexA = (-1); ((indexA < 2) && !hasSameFractionalCoordinate); ++indexA)
							{
								for (int indexB = (-1); ((indexB < 2) && !hasSameFractionalCoordinate); ++indexB)
								{
									for (int indexC = (-1); ((indexC < 2) && !hasSameFractionalCoordinate); ++indexC)
									{
										NumericalVector translatedFractionalCoordinate = originalFractionalCoordinate;
										{
											translatedFractionalCoordinate[0] += static_cast<double>(indexA);
											translatedFractionalCoordinate[1] += static_cast<double>(indexB);
											translatedFractionalCoordinate[2] += static_cast<double>(indexC);
										}


										for (size_type uniqueIndex = 0; uniqueIndex < fractionalCoordinates.size(); ++uniqueIndex)
										{
											if (outputVariables.atomicNumbers.at(index) == atomicNumbers.at(uniqueIndex))
											{
												NumericalVector displacement = outputVariables.unitCell.basisVectors() * (fractionalCoordinates.at(uniqueIndex) - translatedFractionalCoordinate);


												if (displacement.normSquare() < (samePositionDistanceRange * samePositionDistanceRange))
												{
													hasSameFractionalCoordinate = true;
													break;
												}
											}
										}
									}
								}
							}
						}


						if (!hasSameFractionalCoordinate)
						{
							atomicNumbers.push_back(outputVariables.atomicNumbers.at(index));
							fractionalCoordinates.push_back(std::move(originalFractionalCoordinate));
							siteLabels.push_back(outputVariables.siteLabels.at(index));
							wyckoffSymbols.push_back(outputVariables.wyckoffSymbols.at(index));
							siteSymmetrySymbols.push_back(outputVariables.siteSymmetrySymbols.at(index));
							mappingToInputIndices.push_back(outputVariables.mappingToInputIndices.at(index));
						}
					}
				}


				outputVariables.numAtoms = static_cast<size_type>(atomicNumbers.size());
				outputVariables.atomicNumbers = std::move(atomicNumbers);
				outputVariables.fractionalCoordinates = std::move(fractionalCoordinates);
				outputVariables.siteLabels = std::move(siteLabels);
				outputVariables.wyckoffSymbols = std::move(wyckoffSymbols);
				outputVariables.siteSymmetrySymbols = std::move(siteSymmetrySymbols);
				outputVariables.mappingToInputIndices = std::move(mappingToInputIndices);
			}
		}
	}

	return outputVariables;
}

SpaceGroupTypeSearcher::SpglibOutputVariables SpaceGroupTypeSearcher::getConventionalOutputVariables(const SpglibInputVariables& inputVariables) const
{
	SpglibDataset* spglibDatasetPtr = nullptr;
	{
		double prec = _precision;

		for (size_type rep = 0; rep < 20; ++rep)
		{
			spglibDatasetPtr = spg_get_dataset(inputVariables.lattice, inputVariables.positions.get(), inputVariables.types.get(), inputVariables.numAtoms, prec);

			if (spglibDatasetPtr == nullptr)
				prec *= 0.8;
			else
				break;
		}

		if (spglibDatasetPtr == nullptr)
		{
			spg_free_dataset(spglibDatasetPtr);
			throw System::ExceptionServices::ApplicationException{ typeid(*this), "getDataset", "Could not obtain \"spglibDataset\"." };
		}
	}



	try
	{
		const size_type inputNumAtoms = inputVariables.numAtoms;
		SpglibOutputVariables outputVariables;
		{
			outputVariables.spaceGroupNumber = SpaceGroupNumber{ spglibDatasetPtr->spacegroup_number };
			outputVariables.numAtoms = static_cast<size_type>(spglibDatasetPtr->n_std_atoms);
			outputVariables.unitCell.basisVectors() = toNumericalMatrix(spglibDatasetPtr->std_lattice);
		}

		{
			std::map<size_type, size_type> primitiveToInput;
			{
				for (size_type index = 0; index < inputNumAtoms; ++index)
					primitiveToInput.try_emplace(static_cast<size_type>(spglibDatasetPtr->mapping_to_primitive[index]), index);
			}

			for (size_type index = 0; index < outputVariables.numAtoms; ++index)
				outputVariables.mappingToInputIndices.push_back(primitiveToInput.at(static_cast<size_type>(spglibDatasetPtr->std_mapping_to_primitive[index])));
		}


		std::map<size_type, size_type> equivalentAndLabelDictionary;
		{
			std::map<AtomicNumber, std::vector<size_type>> numberAndEquivalentsDictionary;
			{
				for (size_type index = 0; index < outputVariables.numAtoms; ++index)
				{
					size_type equivalentIndex = static_cast<size_type>(spglibDatasetPtr->crystallographic_orbits[outputVariables.mappingToInputIndices.at(index)]);
					AtomicNumber atomicNumber{ spglibDatasetPtr->std_types[index] };


					auto iter = numberAndEquivalentsDictionary.find(atomicNumber);

					if (iter == numberAndEquivalentsDictionary.end())
						numberAndEquivalentsDictionary.emplace(atomicNumber, std::vector<size_type>{ equivalentIndex });
					else
					{
						bool isRegistered = false;
						{
							for (const auto registeredEquivalentIndex : iter->second)
							{
								if (registeredEquivalentIndex == equivalentIndex)
								{
									isRegistered = true;
									break;
								}
							}
						}

						if (!isRegistered)
							iter->second.push_back(equivalentIndex);
					}
				}
			}

			for (const auto& numberAndEquivalents : numberAndEquivalentsDictionary)
			{
				for (size_type index = 0; index < numberAndEquivalents.second.size(); ++index)
					equivalentAndLabelDictionary.emplace(numberAndEquivalents.second.at(index), (1 + index));
			}
		}


		for (int index = 0; index < outputVariables.numAtoms; ++index)
		{
			AtomicNumber atomicNumber{ spglibDatasetPtr->std_types[index] };
			std::string siteLabel = atomicNumber.toElementSymbol();
			siteLabel += std::to_string(equivalentAndLabelDictionary.at(static_cast<size_type>(spglibDatasetPtr->crystallographic_orbits[outputVariables.mappingToInputIndices.at(index)])));

			outputVariables.atomicNumbers.push_back(atomicNumber);
			outputVariables.fractionalCoordinates.push_back(toNumericalVector(spglibDatasetPtr->std_positions[index]));
			outputVariables.siteLabels.push_back(siteLabel);
			outputVariables.wyckoffSymbols.push_back(WyckoffSymbol{ spglibDatasetPtr->wyckoffs[outputVariables.mappingToInputIndices[index]] });
			outputVariables.siteSymmetrySymbols.push_back(SiteSymmetrySymbol{ spglibDatasetPtr->site_symmetry_symbols[outputVariables.mappingToInputIndices[index]] });
		}


		spg_free_dataset(spglibDatasetPtr);
		return outputVariables;
	}


	catch (const System::ExceptionServices::IException& e)
	{
		std::cout << e.toString() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}

	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}
}

SpaceGroupTypeSearcher::SpglibOutputVariables SpaceGroupTypeSearcher::getSymmetrizedInputVariables(const SpglibInputVariables& inputVariables) const
{
	SpglibDataset* spglibDatasetPtr = nullptr;
	{
		double prec = _precision;

		for (size_type rep = 0; rep < 20; ++rep)
		{
			spglibDatasetPtr = spg_get_dataset(inputVariables.lattice, inputVariables.positions.get(), inputVariables.types.get(), inputVariables.numAtoms, prec);

			if (spglibDatasetPtr == nullptr)
				prec *= 0.8;
			else
				break;
		}

		if (spglibDatasetPtr == nullptr)
		{
			spg_free_dataset(spglibDatasetPtr);
			throw System::ExceptionServices::ApplicationException{ typeid(*this), "getDataset", "Could not obtain \"spglibDataset\"." };
		}
	}



	try
	{
		SpglibOutputVariables outputVariables;
		{
			outputVariables.numAtoms = inputVariables.numAtoms;
			outputVariables.unitCell.basisVectors() = toNumericalMatrix(inputVariables.lattice);

			outputVariables.spaceGroupNumber = SpaceGroupNumber{ spglibDatasetPtr->spacegroup_number };
		}


		std::map<size_type, size_type> equivalentAndLabelDictionary;
		{
			std::map<AtomicNumber, std::vector<size_type>> numberAndEquivalentsDictionary;
			{
				for (size_type index = 0; index < outputVariables.numAtoms; ++index)
				{
					size_type equivalentIndex = static_cast<size_type>(spglibDatasetPtr->crystallographic_orbits[index]);
					AtomicNumber atomicNumber{ inputVariables.types[index] };


					auto iter = numberAndEquivalentsDictionary.find(atomicNumber);

					if (iter == numberAndEquivalentsDictionary.end())
						numberAndEquivalentsDictionary.emplace(atomicNumber, std::vector<size_type>{ equivalentIndex });
					else
					{
						bool isRegistered = false;
						{
							for (const auto registeredEquivalentIndex : iter->second)
							{
								if (registeredEquivalentIndex == equivalentIndex)
								{
									isRegistered = true;
									break;
								}
							}
						}

						if (!isRegistered)
							iter->second.push_back(equivalentIndex);
					}
				}
			}

			for (const auto& numberAndEquivalents : numberAndEquivalentsDictionary)
			{
				for (size_type index = 0; index < numberAndEquivalents.second.size(); ++index)
					equivalentAndLabelDictionary.emplace(numberAndEquivalents.second.at(index), (1 + index));
			}
		}


		for (int index = 0; index < outputVariables.numAtoms; ++index)
		{
			AtomicNumber atomicNumber{ inputVariables.types[index] };
			std::string siteLabel = atomicNumber.toElementSymbol();
			siteLabel += std::to_string(equivalentAndLabelDictionary.at(static_cast<size_type>(spglibDatasetPtr->crystallographic_orbits[index])));

			outputVariables.atomicNumbers.push_back(atomicNumber);
			outputVariables.fractionalCoordinates.push_back(toNumericalVector(inputVariables.positions[index]));
			outputVariables.siteLabels.push_back(siteLabel);
			outputVariables.wyckoffSymbols.push_back(WyckoffSymbol{ spglibDatasetPtr->wyckoffs[index] });
			outputVariables.siteSymmetrySymbols.push_back(SiteSymmetrySymbol{ spglibDatasetPtr->site_symmetry_symbols[index] });
		}


		spg_free_dataset(spglibDatasetPtr);
		return outputVariables;
	}


	catch (const System::ExceptionServices::IException& e)
	{
		std::cout << e.toString() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}

	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}
}

std::vector<SpaceGroupTypeSearcher::SymmetryOperation> SpaceGroupTypeSearcher::getSymmetryOperations(const SpglibInputVariables& inputVariables, const SpaceGroupNumber& spaceGroupNumber) const
{
	SpglibDataset* spglibDatasetPtr = nullptr;
	{
		double prec = s_defaultPrecision;

		for (size_type rep = 0; rep < 20; ++rep)
		{
			spglibDatasetPtr = spg_get_dataset(inputVariables.lattice, inputVariables.positions.get(), inputVariables.types.get(), inputVariables.numAtoms, prec);


			if (spglibDatasetPtr == nullptr)
				prec *= 0.8;

			else
			{
				if (SpaceGroupNumber{ spglibDatasetPtr->spacegroup_number } == spaceGroupNumber)
					break;
				else
					prec *= 0.8;
			}
		}

		if (spglibDatasetPtr == nullptr)
		{
			spg_free_dataset(spglibDatasetPtr);
			throw System::ExceptionServices::ApplicationException{ typeid(*this), "getDataset", "Could not obtain \"spglibDataset\"." };
		}
	}


	try
	{
		if (SpaceGroupNumber{ spglibDatasetPtr->spacegroup_number } == spaceGroupNumber)
		{
			std::vector<SymmetryOperation> symmetryOperations;
			{
				for (size_type index = 0; index < spglibDatasetPtr->n_operations; ++index)
					symmetryOperations.push_back(SymmetryOperation{ toNumericalMatrix(spglibDatasetPtr->rotations[index]), toNumericalVector(spglibDatasetPtr->translations[index]) });
			}


			spg_free_dataset(spglibDatasetPtr);
			return symmetryOperations;
		}

		else
			throw System::ExceptionServices::ApplicationException{ typeid(*this), "getDataset", "Could not find the argument space group number." };
	}


	catch (const System::ExceptionServices::IException& e)
	{
		std::cout << e.toString() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}

	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;

		spg_free_dataset(spglibDatasetPtr);
		throw;
	}
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private utility

SpaceGroupTypeSearcher::SpglibInputVariables SpaceGroupTypeSearcher::toSpglibInputVariables(const SpglibOutputVariables& outputVariables) const
{
	SpglibInputVariables spglibInputVariables;
	{
		spglibInputVariables.numAtoms = static_cast<int>(outputVariables.numAtoms);

		spglibInputVariables.lattice[0][0] = outputVariables.unitCell.basisVectors()(0, 0);
		spglibInputVariables.lattice[0][1] = outputVariables.unitCell.basisVectors()(0, 1);
		spglibInputVariables.lattice[0][2] = outputVariables.unitCell.basisVectors()(0, 2);
		spglibInputVariables.lattice[1][0] = outputVariables.unitCell.basisVectors()(1, 0);
		spglibInputVariables.lattice[1][1] = outputVariables.unitCell.basisVectors()(1, 1);
		spglibInputVariables.lattice[1][2] = outputVariables.unitCell.basisVectors()(1, 2);
		spglibInputVariables.lattice[2][0] = outputVariables.unitCell.basisVectors()(2, 0);
		spglibInputVariables.lattice[2][1] = outputVariables.unitCell.basisVectors()(2, 1);
		spglibInputVariables.lattice[2][2] = outputVariables.unitCell.basisVectors()(2, 2);


		spglibInputVariables.types = std::unique_ptr<int[]>{ new int[4 * spglibInputVariables.numAtoms] };
		{
			for (int index = 0; index < spglibInputVariables.numAtoms; ++index)
				spglibInputVariables.types[index] = static_cast<int>(outputVariables.atomicNumbers.at(index));
		}


		spglibInputVariables.positions = std::unique_ptr<double[][3]>{ new double[4 * spglibInputVariables.numAtoms][3] };
		{
			for (int index = 0; index < spglibInputVariables.numAtoms; ++index)
			{
				spglibInputVariables.positions[index][0] = outputVariables.fractionalCoordinates.at(index)[0];
				spglibInputVariables.positions[index][1] = outputVariables.fractionalCoordinates.at(index)[1];
				spglibInputVariables.positions[index][2] = outputVariables.fractionalCoordinates.at(index)[2];
			}
		}
	}

	return spglibInputVariables;
}

// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
