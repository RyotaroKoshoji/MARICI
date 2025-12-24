#ifndef CHEMTOOLKIT_CRYSTALLOGRAPHY_IO_CIFSTREAMWRITER_H
#define CHEMTOOLKIT_CRYSTALLOGRAPHY_IO_CIFSTREAMWRITER_H

#include "FileStream.h"
#include "StreamWriter.h"

#include "CrystallographicStructure.h"
#include "CrystallographicInformation.h"
#include "CrystalStructure.h"


namespace ChemToolkit
{
	namespace Crystallography
	{
		namespace IO
		{
			class CifStreamWriter :private System::IO::FileStream
			{
				using StreamWriter = System::IO::StreamWriter;

				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
				using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;

				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				CifStreamWriter(const std::filesystem::path&, const FileMode&);
				virtual ~CifStreamWriter() = default;

				CifStreamWriter(const CifStreamWriter&) = default;
				CifStreamWriter(CifStreamWriter&&) noexcept = default;
				CifStreamWriter& operator=(const CifStreamWriter&) = default;
				CifStreamWriter& operator=(CifStreamWriter&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Static methods

				static std::string correctExtension() noexcept;
				static bool isCorrectExtension(const std::filesystem::path filePath) noexcept;

			// Static methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void writeCrystallographicInformation(const CrystallographicInformation&) const;

				template <typename A>
				void writeCrystallographicStructure(const CrystallographicStructure<A>&) const;

				template <typename A>
				void writeCrystallographicStructure(const CrystallographicStructure<A>&, const double spaceGroupPrecision) const;

				template <typename A>
				void writeCrystalStructure(const CrystalStructure<A>&) const;

				template <typename A>
				void writeCrystalStructure(const CrystalStructure<A>&, const double spaceGroupPrecision) const;

				template <typename A>
				void writeConventionalCrystalStructure(const CrystalStructure<A>&) const;

				template <typename A>
				void writeConventionalCrystalStructure(const CrystalStructure<A>&, const double spaceGroupPrecision) const;

				template <typename A>
				void writePrimitiveCrystalStructure(const CrystalStructure<A>&) const;

				template <typename A>
				void writePrimitiveCrystalStructure(const CrystalStructure<A>&, const double spaceGroupPrecision) const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void writeCreationDate(StreamWriter&) const;
				void writeChemicalComposition(StreamWriter&, const ChemicalComposition&) const;

				void writeSpaceGroupInformation(StreamWriter&, const CrystallographicInformation&) const;
				void writeUnitCell(StreamWriter&, const CrystallographicInformation&) const;
				void writeAtomSites(StreamWriter&, const CrystallographicInformation&) const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

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
// Static methods

inline std::string ChemToolkit::Crystallography::IO::CifStreamWriter::correctExtension() noexcept
{
	return std::string{ ".cif" };
}

inline bool ChemToolkit::Crystallography::IO::CifStreamWriter::isCorrectExtension(const std::filesystem::path filePath) noexcept
{
	if (filePath.extension().string() == correctExtension())
		return true;
	else
		return false;
}

// Static methods
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

inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeCrystallographicInformation(const CrystallographicInformation& crystallographicInformation) const
{
	StreamWriter streamWriter;
	{
		writeCreationDate(streamWriter);
		writeChemicalComposition(streamWriter, crystallographicInformation.toChemicalComposition());
		writeUnitCell(streamWriter, crystallographicInformation);
		writeSpaceGroupInformation(streamWriter, crystallographicInformation);
		writeAtomSites(streamWriter, crystallographicInformation);
	}

	FileStream::write(streamWriter.allLines());
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeCrystallographicStructure(const CrystallographicStructure<A>& structure) const
{
	writeCrystallographicInformation(structure.toCrystallographicInformation());
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeCrystallographicStructure(const CrystallographicStructure<A>& structure, const double spaceGroupPrecision) const
{
	writeCrystallographicInformation(structure.toCrystallographicInformation(spaceGroupPrecision));
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeCrystalStructure(const CrystalStructure<A>& structure) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.updateSymmetryInformation();

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation());
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeCrystalStructure(const CrystalStructure<A>& structure, const double spaceGroupPrecision) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.updateSymmetryInformation();

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation(spaceGroupPrecision));
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeConventionalCrystalStructure(const CrystalStructure<A>& structure) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.conventionalizeStructure();

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation());
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeConventionalCrystalStructure(const CrystalStructure<A>& structure, const double spaceGroupPrecision) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.conventionalizeStructure(spaceGroupPrecision);

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation(spaceGroupPrecision));
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writePrimitiveCrystalStructure(const CrystalStructure<A>& structure) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.simplifyStructure();

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation());
}

template <typename A>
inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writePrimitiveCrystalStructure(const CrystalStructure<A>& structure, const double spaceGroupPrecision) const
{
	CrystallographicStructure<CrystallographicAtom> crystallographicStructure{ structure.unitCell(), structure.atoms() };
	crystallographicStructure.simplifyStructure(spaceGroupPrecision);

	writeCrystallographicInformation(crystallographicStructure.toCrystallographicInformation(spaceGroupPrecision));
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

inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeChemicalComposition(StreamWriter& streamWriter, const ChemicalComposition& chemicalComposition) const
{
	ChemicalComposition reducedChemicalComposition = chemicalComposition.toReducedChemicalComposition();

	streamWriter.write("_chemical_formula_structural	\'");
	{
		for (const auto& numAndCounting : reducedChemicalComposition)
		{
			streamWriter.write(numAndCounting.first.toElementSymbol());
			streamWriter.write(std::to_string(numAndCounting.second));
			streamWriter.write(" ");
		}

		streamWriter.popBack();
		streamWriter.writeLine("\'");
	}

	streamWriter.write("_chemical_formula_sum	\'");
	{
		for (const auto& numAndCounting : chemicalComposition)
		{
			streamWriter.write(numAndCounting.first.toElementSymbol());
			streamWriter.write(std::to_string(numAndCounting.second));
			streamWriter.write(" ");
		}

		streamWriter.popBack();
		streamWriter.writeLine("\'");
	}
}

inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeSpaceGroupInformation(StreamWriter& streamWriter, const CrystallographicInformation& crystallographicInformation) const
{
	streamWriter.write("_symmetry_Int_Tables_number		");
	streamWriter.writeLine(crystallographicInformation.spaceGroupNumber());

	streamWriter.writeLine("loop_");
	streamWriter.writeLine(" _symmetry_equiv_pos_site_id");
	streamWriter.writeLine(" _symmetry_equiv_pos_as_xyz");


	std::size_t index = 0;

	for (const auto& symmetryOperation : crystallographicInformation.symmetryOperations())
	{
		streamWriter.write(++index);
		streamWriter.write("	");
		streamWriter.writeLine(symmetryOperation.toString());
	}
}

inline void ChemToolkit::Crystallography::IO::CifStreamWriter::writeAtomSites(StreamWriter& streamWriter, const CrystallographicInformation& crystallographicInformation) const
{
	streamWriter.writeLine("loop_");
	streamWriter.writeLine("_atom_site_label");
	streamWriter.writeLine("_atom_site_type_symbol");
	streamWriter.writeLine("_atom_site_symmetry_multiplicity");
	streamWriter.writeLine("_atom_site_Wyckoff_symbol");
	streamWriter.writeLine("_atom_site_fract_x");
	streamWriter.writeLine("_atom_site_fract_y");
	streamWriter.writeLine("_atom_site_fract_z");
	streamWriter.writeLine("_atom_site_occupancy");


	for (const auto& atomSite : crystallographicInformation.atomSites())
	{
		streamWriter.write(atomSite.siteLabel());
		streamWriter.write("  ");
		streamWriter.write(atomSite.toSiteTypeSymbol());
		streamWriter.write("  ");
		streamWriter.write(atomSite.siteMultiplicity());
		streamWriter.write("  ");
		streamWriter.write(atomSite.wyckoffSymbol());
		streamWriter.write("  ");
		streamWriter.write(atomSite.fractionalCoordinate()[0]);
		streamWriter.write("  ");
		streamWriter.write(atomSite.fractionalCoordinate()[1]);
		streamWriter.write("  ");
		streamWriter.write(atomSite.fractionalCoordinate()[2]);
		streamWriter.write("  ");
		streamWriter.writeLine(atomSite.siteOccupancy());
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

#endif // !CHEMTOOLKIT_CRYSTALLOGRAPHY_IO_CIFSTREAMWRITER_H
