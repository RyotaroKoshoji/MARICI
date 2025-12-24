#include "OptimalCrystalStructure.h"

#include <map>
#include <vector>
#include <regex>

#include "InvalidOperationException.h"

#include "InvalidFileException.h"
#include "StreamWriter.h"

#include "ConstrainingCrystalStructure.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

OptimalCrystalStructure::OptimalCrystalStructure() noexcept
	: CrystallographicStructure{}
{
}

OptimalCrystalStructure::OptimalCrystalStructure(const UnitCell& cell, const std::vector<OptimalAtom>& atoms) noexcept
	: CrystallographicStructure{ cell, atoms }
{
}

OptimalCrystalStructure::OptimalCrystalStructure(const UnitCell& cell, std::vector<OptimalAtom>&& atoms) noexcept
	: CrystallographicStructure{ cell, std::move(atoms) }
{
}

OptimalCrystalStructure::OptimalCrystalStructure(const CrystallographicStructure<ChemToolkit::Crystallography::CrystallographicAtom>& structure)
	: CrystallographicStructure{ structure.unitCell() }
{
	std::vector<OptimalAtom> atoms;
	{
		for (const auto& atom : structure.atoms())
			atoms.push_back(OptimalAtom{ atom.ionicAtomicNumber(), atom.cartesianCoordinate() });
	}

	CrystallographicStructure::setAtoms(std::move(atoms));
}

OptimalCrystalStructure::OptimalCrystalStructure(const ConstrainingCrystalStructure& structure)
	: CrystallographicStructure{ structure.unitCell() }
{
	std::vector<OptimalAtom> atomicArrangement;
	{
		for (const auto& constrainingAtom : structure.atoms())
			atomicArrangement.push_back(OptimalAtom{ constrainingAtom });


		for (size_type index = 0; index < atomicArrangement.size(); ++index)
		{
			for (const auto& originalAtomIndex : structure.atoms()[index].getCovalentBondedOriginalAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(originalAtomIndex).ionicAtomicNumber());

			for (const auto& originalAtomIndex : structure.atoms()[index].getIonicBondedOriginalAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(originalAtomIndex).ionicAtomicNumber());

			for (const auto& translatedAtomIndex : structure.atoms()[index].getCovalentBondedTranslatedAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(translatedAtomIndex.originalIndex()).ionicAtomicNumber());

			for (const auto& translatedAtomIndex : structure.atoms()[index].getIonicBondedTranslatedAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(translatedAtomIndex.originalIndex()).ionicAtomicNumber());
		}


		for (const auto& coordinationPolyhedraLinking : structure.getCoordinationPolyhedraLinkings())
		{
			if (1 == coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				const OptimalAtom& vertexAtom = atomicArrangement.at(coordinationPolyhedraLinking.second.back().originalIndex());
				IonicAtomicNumber vertexAtomicNumber{ vertexAtom.ionicAtomicNumber() };

				originalAtom.addVertexSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), vertexAtomicNumber));
				translatedAtom.addVertexSharings(std::make_pair(originalAtom.ionicAtomicNumber(), vertexAtomicNumber));
			}

			else if (2 == coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				std::array<IonicAtomicNumber, 2> edgeAtomicNumbers;
				{
					edgeAtomicNumbers[0] = atomicArrangement[coordinationPolyhedraLinking.second[0].originalIndex()].ionicAtomicNumber();
					edgeAtomicNumbers[1] = atomicArrangement[coordinationPolyhedraLinking.second[1].originalIndex()].ionicAtomicNumber();
				}

				originalAtom.addEdgeSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), edgeAtomicNumbers));
				translatedAtom.addEdgeSharings(std::make_pair(originalAtom.ionicAtomicNumber(), edgeAtomicNumbers));
			}

			else if (2 < coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				std::vector<IonicAtomicNumber> faceAtomicNumbers;
				{
					for (const auto& translatedIndex : coordinationPolyhedraLinking.second)
						faceAtomicNumbers.push_back(atomicArrangement[translatedIndex.originalIndex()].ionicAtomicNumber());
				}

				originalAtom.addFaceSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), faceAtomicNumbers));
				translatedAtom.addFaceSharings(std::make_pair(originalAtom.ionicAtomicNumber(), faceAtomicNumbers));
			}

			else
				throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "import", "The number of common bridging atoms is zero." };
		}
	}

	CrystallographicStructure::setAtoms(std::move(atomicArrangement));
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
// Methods

void OptimalCrystalStructure::initialize() noexcept
{
	CrystallographicStructure::initialize();
}

void OptimalCrystalStructure::import(const ConstrainingCrystalStructure& structure)
{
	initialize();
	
	std::vector<OptimalAtom> atomicArrangement;
	{
		for (const auto& constrainingAtom : structure.atoms())
			atomicArrangement.push_back(OptimalAtom{ constrainingAtom });


		for (size_type index = 0; index < atomicArrangement.size(); ++index)
		{
			for (const auto& originalAtomIndex : structure.atoms()[index].getCovalentBondedOriginalAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(originalAtomIndex).ionicAtomicNumber());

			for (const auto& originalAtomIndex : structure.atoms()[index].getIonicBondedOriginalAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(originalAtomIndex).ionicAtomicNumber());

			for (const auto& translatedAtomIndex : structure.atoms()[index].getCovalentBondedTranslatedAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(translatedAtomIndex.originalIndex()).ionicAtomicNumber());

			for (const auto& translatedAtomIndex : structure.atoms()[index].getIonicBondedTranslatedAtomIndices())
				atomicArrangement[index].addCoordination(structure.atoms().at(translatedAtomIndex.originalIndex()).ionicAtomicNumber());
		}


		for (const auto& coordinationPolyhedraLinking : structure.getCoordinationPolyhedraLinkings())
		{
			if (1 == coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				const OptimalAtom& vertexAtom = atomicArrangement.at(coordinationPolyhedraLinking.second.back().originalIndex());
				IonicAtomicNumber vertexAtomicNumber{ vertexAtom.ionicAtomicNumber() };

				originalAtom.addVertexSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), vertexAtomicNumber));
				translatedAtom.addVertexSharings(std::make_pair(originalAtom.ionicAtomicNumber(), vertexAtomicNumber));
			}

			else if (2 == coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				std::array<IonicAtomicNumber, 2> edgeAtomicNumbers;
				{
					edgeAtomicNumbers[0] = atomicArrangement[coordinationPolyhedraLinking.second[0].originalIndex()].ionicAtomicNumber();
					edgeAtomicNumbers[1] = atomicArrangement[coordinationPolyhedraLinking.second[1].originalIndex()].ionicAtomicNumber();
				}

				originalAtom.addEdgeSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), edgeAtomicNumbers));
				translatedAtom.addEdgeSharings(std::make_pair(originalAtom.ionicAtomicNumber(), edgeAtomicNumbers));
			}

			else if (2 < coordinationPolyhedraLinking.second.size())
			{
				OptimalAtom& originalAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.originalAtomIndex());
				OptimalAtom& translatedAtom = atomicArrangement.at(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex());

				std::vector<IonicAtomicNumber> faceAtomicNumbers;
				{
					for (const auto& translatedIndex : coordinationPolyhedraLinking.second)
						faceAtomicNumbers.push_back(atomicArrangement[translatedIndex.originalIndex()].ionicAtomicNumber());
				}

				originalAtom.addFaceSharings(std::make_pair(translatedAtom.ionicAtomicNumber(), faceAtomicNumbers));
				translatedAtom.addFaceSharings(std::make_pair(originalAtom.ionicAtomicNumber(), faceAtomicNumbers));
			}

			else
				throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "import", "The number of common bridging atoms is zero." };
		}
	}

	setUnitCell(ChemToolkit::Crystallography::UnitCell{ structure.unitCell().basisVectors() });
	setAtoms(std::move(atomicArrangement));
}

std::string OptimalCrystalStructure::toAbridgedStructuralFingerprint() const
{
	std::vector<OptimalAtom> atomicArrangement{ atoms() };
	std::sort(atomicArrangement.begin(), atomicArrangement.end());

	std::map<AtomicNumber, size_type> numberAndCount;
	System::IO::StreamWriter streamWriter;
	{
		for (const auto& atom : atomicArrangement)
		{
			AtomicNumber atomicNumber = atom.ionicAtomicNumber().atomicNumber();
			{
				auto iter = numberAndCount.find(atomicNumber);

				if (iter == numberAndCount.end())
					numberAndCount.emplace(atomicNumber, 1);
				else
					iter->second += 1;
			}

			streamWriter.write("$");
			streamWriter.write(atomicNumber.toElementSymbol());
			streamWriter.write(numberAndCount[atomicNumber]);
			streamWriter.writeLine("_environment");

			streamWriter.writeLine(atom.getAbridgedFingerprint());
			streamWriter.breakLine();
		}
	}

	return streamWriter.allTexts();
}

std::string OptimalCrystalStructure::toStructuralFingerprint() const
{
	System::IO::StreamWriter streamWriter;
	{
		std::vector<std::pair<OptimalAtom, size_type>> uniqueOptimalAtoms = getUniqueOptimalAtoms();
		std::map<AtomicNumber, size_type> numberAndCount;


		streamWriter.write("Space.Group.Number\t");
		streamWriter.writeLine(spaceGroupNumber());
		streamWriter.breakLine();

		for (const auto& uniqueOptimalAtom : uniqueOptimalAtoms)
		{
			AtomicNumber atomicNumber = uniqueOptimalAtom.first.ionicAtomicNumber().atomicNumber();
			{
				auto iter = numberAndCount.find(atomicNumber);

				if (iter == numberAndCount.end())
					numberAndCount.emplace(atomicNumber, 1);
				else
					iter->second += 1;
			}


			streamWriter.write("&");
			streamWriter.write(atomicNumber.toElementSymbol());
			streamWriter.write(numberAndCount[atomicNumber]);
			streamWriter.writeLine("_environment");

			streamWriter.writeLine(uniqueOptimalAtom.first.getFingerprint());
			streamWriter.write("\tSite.Multiplicity\t\t");
			streamWriter.writeLine(uniqueOptimalAtom.second);
			streamWriter.breakLine();
		}
	}

	return streamWriter.allTexts();
}

std::string OptimalCrystalStructure::toDetailedStructuralFingerprint() const
{
	System::IO::StreamWriter streamWriter;
	{
		std::vector<std::pair<OptimalAtom, size_type>> uniqueOptimalAtoms = getUniqueOptimalAtoms();
		std::map<AtomicNumber, size_type> numberAndCount;


		streamWriter.write("Space.Group.Number\t");
		streamWriter.writeLine(spaceGroupNumber());
		streamWriter.breakLine();

		for (const auto& uniqueOptimalAtom : uniqueOptimalAtoms)
		{
			AtomicNumber atomicNumber = uniqueOptimalAtom.first.ionicAtomicNumber().atomicNumber();
			{
				auto iter = numberAndCount.find(atomicNumber);

				if (iter == numberAndCount.end())
					numberAndCount.emplace(atomicNumber, 1);
				else
					iter->second += 1;
			}


			streamWriter.write("&");
			streamWriter.write(atomicNumber.toElementSymbol());
			streamWriter.write(numberAndCount[atomicNumber]);
			streamWriter.writeLine("_environment");

			streamWriter.writeLine(uniqueOptimalAtom.first.getDetailedFingerprint());
			streamWriter.write("\tSite.Multiplicity\t\t");
			streamWriter.writeLine(uniqueOptimalAtom.second);
			streamWriter.breakLine();
		}
	}

	return streamWriter.allTexts();
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

std::vector<std::pair<OptimalAtom, OptimalCrystalStructure::size_type>> OptimalCrystalStructure::getUniqueOptimalAtoms() const
{
	std::vector<std::pair<OptimalAtom, size_type>> uniqueOptimalAtoms;
	{
		std::vector<OptimalAtom> atomicArrangement = atoms();
		std::sort(atomicArrangement.begin(), atomicArrangement.end());


		for (const auto& atom : atomicArrangement)
		{
			auto uniqueIter = uniqueOptimalAtoms.begin();
			{
				for (; uniqueIter != uniqueOptimalAtoms.end(); ++uniqueIter)
				{
					if (uniqueIter->first == atom)
						break;
				}
			}

			if (uniqueIter == uniqueOptimalAtoms.end())
				uniqueOptimalAtoms.push_back(std::make_pair(atom, 1));
			else
				uniqueIter->second += 1;
		}
	}

	return uniqueOptimalAtoms;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
