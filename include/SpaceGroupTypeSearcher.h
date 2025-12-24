#ifndef CHEMTOOLKIT_CRYSTALLOGRAPHY_SYMMETRY_SPACEGROUPTYPESEARCHER_H
#define CHEMTOOLKIT_CRYSTALLOGRAPHY_SYMMETRY_SPACEGROUPTYPESEARCHER_H

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "ArgumentOutOfRangeException.h"

#include "NumericalVector.h"
#include "NumericalMatrix.h"

#include "AtomicNumber.h"
#include "CrystallographicSymmetryOperation.h"
#include "SpaceGroupNumber.h"
#include "SiteSymmetrySymbol.h"
#include "WyckoffSymbol.h"

#include "UnitCell.h"
#include "CrystallographicInformation.h"


namespace ChemToolkit
{
	namespace Crystallography
	{
		namespace Symmetry
		{
			class SpaceGroupTypeSearcher
			{
				using size_type = unsigned short;
				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
				using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;
				using SymmetryOperation = ChemToolkit::Crystallography::Symmetry::CrystallographicSymmetryOperation;

				struct SpglibInputVariables
				{
					int numAtoms{ 0 };
					double lattice[3][3];
					std::unique_ptr<int[]> types;
					std::unique_ptr<double[][3]> positions;
					std::vector<size_type> mappingToInputIndices;
				};

				struct SpglibOutputVariables
				{
					SpaceGroupNumber spaceGroupNumber;
					size_type numAtoms{ 0 };
					UnitCell unitCell;

					std::vector<AtomicNumber> atomicNumbers;
					std::vector<NumericalVector> fractionalCoordinates;

					std::vector<std::string> siteLabels;
					std::vector<WyckoffSymbol> wyckoffSymbols;
					std::vector<SiteSymmetrySymbol> siteSymmetrySymbols;
					std::vector<size_type> mappingToInputIndices;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				SpaceGroupTypeSearcher() noexcept;
				explicit SpaceGroupTypeSearcher(const double precision);

				virtual ~SpaceGroupTypeSearcher() = default;

				SpaceGroupTypeSearcher(const SpaceGroupTypeSearcher&) = default;
				SpaceGroupTypeSearcher(SpaceGroupTypeSearcher&&) noexcept = default;
				SpaceGroupTypeSearcher& operator=(const SpaceGroupTypeSearcher&) = default;
				SpaceGroupTypeSearcher& operator=(SpaceGroupTypeSearcher&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				double precision() const noexcept;
				static double defaultPrecision() noexcept;

				void setPrecision() noexcept;
				void setPrecision(const double);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				template <typename A>
				void applyDelaunayReduction(UnitCell&, std::vector<A>&) const;

				template <typename A>
				void toPrimitive(UnitCell&, std::vector<A>&, SpaceGroupNumber&) const;

				template <typename A>
				void conventionalize(UnitCell&, std::vector<A>&, SpaceGroupNumber&) const;

				template <typename A>
				void updateSymmetryInformation(const UnitCell&, std::vector<A>&, SpaceGroupNumber&) const;

				template <typename A>
				CrystallographicInformation getCrystallographicInformation(const UnitCell&, const std::vector<A>&, const SpaceGroupNumber) const;

				template <typename A>
				CrystallographicInformation getCrystallographicInformation(const UnitCell&, const std::vector<A>&) const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				UnitCell getDelaunayReducedUnitCell(const UnitCell&) const;
				SpglibOutputVariables getPrimitiveOutputVariables(const SpglibInputVariables&) const;
				SpglibOutputVariables getConventionalOutputVariables(const SpglibInputVariables&) const;
				SpglibOutputVariables getSymmetrizedInputVariables(const SpglibInputVariables&) const;
				std::vector<SymmetryOperation> getSymmetryOperations(const SpglibInputVariables&, const SpaceGroupNumber&) const;


				template <typename A>
				SpglibInputVariables toSpglibInputVariables(const UnitCell&, const std::vector<A>&) const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private utility

				NumericalVector toNumericalVector(const double vec[3]) const;
				NumericalMatrix toNumericalMatrix(const double mat[3][3]) const;
				NumericalMatrix toNumericalMatrix(const int mat[3][3]) const;

				SpglibInputVariables toSpglibInputVariables(const SpglibOutputVariables&) const;
				void removeEquivalentAtoms(SpglibOutputVariables&, const UnitCell& prevUnitCell) const;

			// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				double _precision;

				static double s_defaultPrecision;
			};
		}
	}
}

// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Property

inline double ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::precision() const noexcept
{
	return _precision;
}

inline double ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::defaultPrecision() noexcept
{
	return s_defaultPrecision;
}

inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::setPrecision() noexcept
{
	_precision = s_defaultPrecision;
}

inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::setPrecision(const double val)
{
	if (0.0 < val)
		_precision = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), " setPrecision", "Input precision is not more than zero." };
}

// Property
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
// Methods

template <typename A>
inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::applyDelaunayReduction(UnitCell& unitCell, std::vector<A>& atomicArrangement) const
{
	unitCell = getDelaunayReducedUnitCell(unitCell);
	NumericalMatrix inverseBasisVectors = unitCell.getInverseBasisVectors();

	for (auto& atom : atomicArrangement)
	{
		NumericalVector fractionalCoordinate = inverseBasisVectors * atom.cartesianCoordinate();
		{
			int floorValueA = static_cast<int>(std::floor(fractionalCoordinate[0]));
			int floorValueB = static_cast<int>(std::floor(fractionalCoordinate[1]));
			int floorValueC = static_cast<int>(std::floor(fractionalCoordinate[2]));

			if (floorValueA != 0)
				fractionalCoordinate[0] -= static_cast<double>(floorValueA);
			if (floorValueB != 0)
				fractionalCoordinate[1] -= static_cast<double>(floorValueB);
			if (floorValueC != 0)
				fractionalCoordinate[2] -= static_cast<double>(floorValueC);
		}

		atom.cartesianCoordinate() = unitCell.basisVectors() * fractionalCoordinate;
	}
}

template <typename A>
inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::toPrimitive(UnitCell& unitCell, std::vector<A>& atoms, SpaceGroupNumber& spaceGroupNumber) const
{
	SpglibOutputVariables variables = getPrimitiveOutputVariables(toSpglibInputVariables(unitCell, atoms));

	using Atom = A;
	std::vector<Atom> oldAtoms = std::move(atoms);
	atoms.clear();
	{
		for (size_type index = 0; index < static_cast<size_type>(variables.numAtoms); ++index)
		{
			Atom atom = oldAtoms.at(variables.mappingToInputIndices.at(index));
			atom.cartesianCoordinate() = (variables.unitCell.basisVectors() * variables.fractionalCoordinates.at(index));
			atom.setSiteLabel(variables.siteLabels.at(index));
			atom.setSiteSymmetrySymbol(variables.siteSymmetrySymbols.at(index));
			atom.setWyckoffSymbol(variables.wyckoffSymbols.at(index));

			atoms.push_back(std::move(atom));
		}
	}


	unitCell = variables.unitCell;
	spaceGroupNumber = variables.spaceGroupNumber;
}

template <typename A>
inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::conventionalize(UnitCell& unitCell, std::vector<A>& atoms, SpaceGroupNumber& spaceGroupNumber) const
{
	SpglibOutputVariables variables = getConventionalOutputVariables(toSpglibInputVariables(unitCell, atoms));

	using Atom = A;
	std::vector<Atom> oldAtoms = std::move(atoms);
	atoms.clear();
	{
		for (size_type index = 0; index < static_cast<size_type>(variables.numAtoms); ++index)
		{
			Atom atom = oldAtoms.at(variables.mappingToInputIndices.at(index));
			atom.cartesianCoordinate() = (variables.unitCell.basisVectors() * variables.fractionalCoordinates.at(index));
			atom.setSiteLabel(variables.siteLabels.at(index));
			atom.setSiteSymmetrySymbol(variables.siteSymmetrySymbols.at(index));
			atom.setWyckoffSymbol(variables.wyckoffSymbols.at(index));

			atoms.push_back(std::move(atom));
		}
	}


	unitCell = variables.unitCell;
	spaceGroupNumber = variables.spaceGroupNumber;
}

template <typename A>
inline void ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::updateSymmetryInformation(const UnitCell& unitCell, std::vector<A>& atoms, SpaceGroupNumber& spaceGroupNumber) const
{
	SpglibOutputVariables variables = getSymmetrizedInputVariables(toSpglibInputVariables(unitCell, atoms));


	spaceGroupNumber = variables.spaceGroupNumber;

	for (size_type index = 0; index < atoms.size(); ++index)
	{
		atoms[index].setSiteLabel(variables.siteLabels.at(index));
		atoms[index].setSiteSymmetrySymbol(variables.siteSymmetrySymbols.at(index));
		atoms[index].setWyckoffSymbol(variables.wyckoffSymbols.at(index));
	}
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicInformation ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::getCrystallographicInformation(const UnitCell& unitCell, const std::vector<A>& atoms, const SpaceGroupNumber spaceGroupNumber) const
{
	std::vector<SymmetryOperation> symmetryOperations = getSymmetryOperations(toSpglibInputVariables(unitCell, atoms), spaceGroupNumber);
	std::vector<AtomSite> atomSites;
	{
		const NumericalMatrix inverseBasisVectors = unitCell.getInverseBasisVectors();


		for (const auto& atom : atoms)
		{
			if (atom.hasSymmetryInformation())
			{
				auto siteIter = atomSites.begin();
				{
					for (; siteIter != atomSites.end(); ++siteIter)
					{
						if (siteIter->siteLabel() == atom.siteLabel())
						{
							size_type multiplicity = siteIter->siteMultiplicity();
							siteIter->setSiteMultiplicity(++multiplicity);
							break;
						}
					}
				}

				if (siteIter == atomSites.end())
				{
					AtomSite atomSite;
					{
						atomSite.setAtomicNumber(atom.ionicAtomicNumber().atomicNumber());
						atomSite.setFormalCharge(atom.ionicAtomicNumber().formalCharge());
						atomSite.setFractionalCoordinate((inverseBasisVectors * atom.cartesianCoordinate()));

						atomSite.setSiteLabel(atom.siteLabel());
						atomSite.setSiteSymmetrySymbol(atom.siteSymmetrySymbol());
						atomSite.setWyckoffSymbol(atom.wyckoffSymbol());

						atomSite.setSiteOccupancy(1.0);
						atomSite.setSiteMultiplicity(1);
					}

					atomSites.push_back(std::move(atomSite));
				}
			}

			else
				throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "getCrystallographicInformation", "Argument atoms do not have symmetry information." };
		}
	}


	CrystallographicInformation crystallographicInformation;
	crystallographicInformation.setUnitCell(unitCell);
	crystallographicInformation.setAtomSites(std::move(atomSites));
	crystallographicInformation.setSpaceGroupNumber(spaceGroupNumber);
	crystallographicInformation.setSymmetryOperations(std::move(symmetryOperations));

	return crystallographicInformation;
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicInformation ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::getCrystallographicInformation(const UnitCell& unitCell, const std::vector<A>& atoms) const
{
	std::vector<SymmetryOperation> symmetryOperations;
	symmetryOperations.push_back(SymmetryOperation{ MathToolkit::LinearAlgebra::getIdentityMatrix<double, 3>() });


	std::vector<AtomSite> atomSites;
	{
		const NumericalMatrix inverseBasisVectors = unitCell.getInverseBasisVectors();
		std::map<size_type, size_type> siteLabelDictionary;
		{
			std::map<AtomicNumber, size_type> numberAndCount;
			{
				for (size_type index = 0; index < atoms.size(); ++index)
				{
					AtomicNumber atomicNumber = atoms[index].ionicAtomicNumber().atomicNumber();
					auto iter = numberAndCount.find(atomicNumber);


					if (iter == numberAndCount.end())
					{
						numberAndCount.emplace(atomicNumber, 1);
						siteLabelDictionary.emplace(index, 1);
					}

					else
					{
						iter->second += 1;
						siteLabelDictionary.emplace(index, iter->second);
					}
				}
			}
		}


		for (size_type index = 0; index < atoms.size(); ++index)
		{
			std::string siteLabel = atoms[index].ionicAtomicNumber().atomicNumber().toElementSymbol();
			siteLabel += std::to_string(siteLabelDictionary.at(index));


			AtomSite atomSite;
			{
				atomSite.setAtomicNumber(atoms[index].ionicAtomicNumber().atomicNumber());
				atomSite.setFormalCharge(atoms[index].ionicAtomicNumber().formalCharge());
				atomSite.setFractionalCoordinate((inverseBasisVectors * atoms[index].cartesianCoordinate()));

				atomSite.setSiteLabel(siteLabel);
				atomSite.setSiteSymmetrySymbol(SiteSymmetrySymbol{ "1" });
				atomSite.setWyckoffSymbol(WyckoffSymbol{ "a" });

				atomSite.setSiteOccupancy(1.0);
				atomSite.setSiteMultiplicity(1);
			}

			atomSites.push_back(std::move(atomSite));
		}
	}


	CrystallographicInformation crystallographicInformation;
	crystallographicInformation.setUnitCell(unitCell);
	crystallographicInformation.setAtomSites(std::move(atomSites));
	crystallographicInformation.setSpaceGroupNumber(SpaceGroupNumber{ 1 });
	crystallographicInformation.setSymmetryOperations(std::move(symmetryOperations));

	return crystallographicInformation;
}

// Methods
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

template <typename A>
inline ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::SpglibInputVariables ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::toSpglibInputVariables(const UnitCell& cell, const std::vector<A>& atoms) const
{
	SpglibInputVariables spglibInputVariables;
	{
		spglibInputVariables.numAtoms = static_cast<int>(atoms.size());

		spglibInputVariables.lattice[0][0] = cell.basisVectors()(0, 0);
		spglibInputVariables.lattice[0][1] = cell.basisVectors()(0, 1);
		spglibInputVariables.lattice[0][2] = cell.basisVectors()(0, 2);
		spglibInputVariables.lattice[1][0] = cell.basisVectors()(1, 0);
		spglibInputVariables.lattice[1][1] = cell.basisVectors()(1, 1);
		spglibInputVariables.lattice[1][2] = cell.basisVectors()(1, 2);
		spglibInputVariables.lattice[2][0] = cell.basisVectors()(2, 0);
		spglibInputVariables.lattice[2][1] = cell.basisVectors()(2, 1);
		spglibInputVariables.lattice[2][2] = cell.basisVectors()(2, 2);

		spglibInputVariables.types = std::unique_ptr<int[]>{ new int[4 * spglibInputVariables.numAtoms] };
		{
			for (int index = 0; index < spglibInputVariables.numAtoms; ++index)
				spglibInputVariables.types[index] = static_cast<int>(atoms.at(index).ionicAtomicNumber().atomicNumber());
		}

		spglibInputVariables.positions = std::unique_ptr<double[][3]>{ new double[4 * spglibInputVariables.numAtoms][3] };
		{
			NumericalMatrix inverseBasisVectors = cell.getInverseBasisVectors();

			for (int index = 0; index < spglibInputVariables.numAtoms; ++index)
			{
				NumericalVector fractionalCoordinate = inverseBasisVectors * atoms.at(index).cartesianCoordinate();
				spglibInputVariables.positions[index][0] = fractionalCoordinate[0];
				spglibInputVariables.positions[index][1] = fractionalCoordinate[1];
				spglibInputVariables.positions[index][2] = fractionalCoordinate[2];
			}
		}
	}

	return spglibInputVariables;
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

inline ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::NumericalVector ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::toNumericalVector(const double rowVec[3]) const
{
	NumericalVector vec;
	{
		vec[0] = rowVec[0];
		vec[1] = rowVec[1];
		vec[2] = rowVec[2];
	}

	return vec;
}

inline ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::NumericalMatrix ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::toNumericalMatrix(const double rowMat[3][3]) const
{
	NumericalMatrix mat;
	{
		mat(0, 0) = rowMat[0][0];
		mat(0, 1) = rowMat[0][1];
		mat(0, 2) = rowMat[0][2];
		mat(1, 0) = rowMat[1][0];
		mat(1, 1) = rowMat[1][1];
		mat(1, 2) = rowMat[1][2];
		mat(2, 0) = rowMat[2][0];
		mat(2, 1) = rowMat[2][1];
		mat(2, 2) = rowMat[2][2];
	}

	return mat;
}

inline ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::NumericalMatrix ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher::toNumericalMatrix(const int rowMat[3][3]) const
{
	NumericalMatrix mat;
	{
		mat(0, 0) = static_cast<double>(rowMat[0][0]);
		mat(0, 1) = static_cast<double>(rowMat[0][1]);
		mat(0, 2) = static_cast<double>(rowMat[0][2]);
		mat(1, 0) = static_cast<double>(rowMat[1][0]);
		mat(1, 1) = static_cast<double>(rowMat[1][1]);
		mat(1, 2) = static_cast<double>(rowMat[1][2]);
		mat(2, 0) = static_cast<double>(rowMat[2][0]);
		mat(2, 1) = static_cast<double>(rowMat[2][1]);
		mat(2, 2) = static_cast<double>(rowMat[2][2]);
	}

	return mat;
}

// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !CHEMTOOLKIT_CRYSTALLOGRAPHY_SYMMETRY_SPACEGROUPTYPESEARCHER_H
