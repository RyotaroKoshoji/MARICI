#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_FEASIBLEPOLYHEDRACONNECTIONSDICTIONARY_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_FEASIBLEPOLYHEDRACONNECTIONSDICTIONARY_H

#include <unordered_map>
#include <utility>

#include "StreamReader.h"

#include "AtomicNumber.h"
#include "ChemicalComposition.h"

#include "IonicAtomicNumber.h"

#include "ConstrainingAtomicSpecies.h"

#include "CoordinationPolyhedraConnectionParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			class FeasiblePolyhedraConnectionsDictionary
			{
				using size_type = unsigned short;
				using CoordinationPolyhedraConnectionParameters = Internal::CoordinationPolyhedraConnectionParameters;
				
				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
				using BasicChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;
				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<IonicAtomicNumber>;

				using ConstrainingAtomicSpecies = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies;
				using ConstrainingChemicalComposition = ChemToolkit::Generic::ChemicalComposition<ConstrainingAtomicSpecies>;


				struct KeyHasher
				{
					std::size_t operator()(const std::pair<IonicAtomicNumber, IonicAtomicNumber>&) const noexcept;
				};

				struct KeyEqual
				{
					bool operator()(const std::pair<IonicAtomicNumber, IonicAtomicNumber>&, const std::pair<IonicAtomicNumber, IonicAtomicNumber>&) const noexcept;
				};

				using FeasibleLinkingDictionary = std::unordered_map<std::pair<IonicAtomicNumber, IonicAtomicNumber>, std::vector<BasicChemicalComposition>, KeyHasher, KeyEqual>;


// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			private:
				FeasiblePolyhedraConnectionsDictionary() noexcept = delete;

			public:
				virtual ~FeasiblePolyhedraConnectionsDictionary() = default;

				FeasiblePolyhedraConnectionsDictionary(const FeasiblePolyhedraConnectionsDictionary&) = default;
				FeasiblePolyhedraConnectionsDictionary(FeasiblePolyhedraConnectionsDictionary&&) noexcept = default;
				FeasiblePolyhedraConnectionsDictionary& operator=(const FeasiblePolyhedraConnectionsDictionary&) = default;
				FeasiblePolyhedraConnectionsDictionary& operator=(FeasiblePolyhedraConnectionsDictionary&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				static void initialize() noexcept;
				static void initialize(const System::IO::StreamReader&);
				static void initialize(const ChemicalComposition& crystalComposition);
				static void initialize(const ConstrainingChemicalComposition& constrainingCrystalComposition);

				static FeasibleLinkingDictionary& getDictionary() noexcept;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				static std::pair<IonicAtomicNumber, IonicAtomicNumber> getIonicAtomicNumberKey(const IonicAtomicNumber&, const IonicAtomicNumber&);

				static ConstrainingChemicalComposition toConstrainingChemicalComposition(const ChemicalComposition&);

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				static CoordinationPolyhedraConnectionParameters s_coordinationPolyhedraConnectionParameters;

				static FeasibleLinkingDictionary s_feasiblePolyhedraLinkingDictionary;
			};



			inline std::size_t FeasiblePolyhedraConnectionsDictionary::KeyHasher::operator()(const std::pair<IonicAtomicNumber, IonicAtomicNumber>& ianp) const noexcept
			{
				return (ianp.first.getHashCode() ^ ianp.second.getHashCode());
			}

			inline bool FeasiblePolyhedraConnectionsDictionary::KeyEqual::operator()(const std::pair<IonicAtomicNumber, IonicAtomicNumber>& ianpa, const std::pair<IonicAtomicNumber, IonicAtomicNumber>& ianpb) const noexcept
			{
				return (ianpa == ianpb);
			}
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
// Methods

inline std::unordered_map<std::pair<MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber, MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber>, std::vector<MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::BasicChemicalComposition>, MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::KeyHasher, MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::KeyEqual>& MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::getDictionary() noexcept
{
	return s_feasiblePolyhedraLinkingDictionary;
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

inline std::pair<MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber, MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber> MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::getIonicAtomicNumberKey(const IonicAtomicNumber& iana, const IonicAtomicNumber& ianb)
{
	if (ianb < iana)
		return std::make_pair(ianb, iana);
	else
		return std::make_pair(iana, ianb);
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::ConstrainingChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::toConstrainingChemicalComposition(const ChemicalComposition& chemicalComposition)
{
	ConstrainingChemicalComposition constrainingChemicalComposition;
	{
		for (const auto& numAndCount : chemicalComposition)
			constrainingChemicalComposition.add(ConstrainingAtomicSpecies{ numAndCount.first }, numAndCount.second);
	}

	return constrainingChemicalComposition;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_FEASIBLEPOLYHEDRACONNECTIONSDICTIONARY_H
