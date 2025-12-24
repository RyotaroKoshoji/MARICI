#ifndef DFTTOOLKIT_OPENMX_IO_OPENMXMDSTREAMWRITER_H
#define DFTTOOLKIT_OPENMX_IO_OPENMXMDSTREAMWRITER_H

#include <string>

#include "NumericalVector.h"
#include "NumericalMatrix.h"

#include "FileStream.h"
#include "StreamWriter.h"

#include "UnitCell.h"
#include "CrystalStructure.h"
#include "CrystallographicStructure.h"


namespace DftToolkit
{
	namespace Openmx
	{
		namespace IO
		{
			class OpenmxMdStreamWriter :private System::IO::FileStream
			{
				using StreamWriter = System::IO::StreamWriter;

				using size_type = std::size_t;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
				using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;

				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using UnitCell = ChemToolkit::Crystallography::UnitCell;

				template <typename A>
				using CrystallographicStructure = ChemToolkit::Crystallography::CrystallographicStructure<A>;

				template <typename A>
				using CrystalStructure = ChemToolkit::Crystallography::CrystalStructure<A>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				OpenmxMdStreamWriter() noexcept;
				OpenmxMdStreamWriter(const std::filesystem::path&, const FileMode&);
				virtual ~OpenmxMdStreamWriter() = default;

				OpenmxMdStreamWriter(const OpenmxMdStreamWriter&) = default;
				OpenmxMdStreamWriter(OpenmxMdStreamWriter&&) noexcept = default;
				OpenmxMdStreamWriter& operator=(const OpenmxMdStreamWriter&) = default;
				OpenmxMdStreamWriter& operator=(OpenmxMdStreamWriter&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Static methods

				static std::string correctExtension() noexcept;
				static bool isCorrectExtension(const std::filesystem::path& filePath);

			// Static methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				template <typename A>
				void record(const CrystallographicStructure<A>&, const size_type elapsedCount);

				template <typename A>
				void record(const CrystallographicStructure<A>&, const double elapsedTime);

				template <typename A>
				void record(const CrystalStructure<A>&, const size_type elapsedCount);

				template <typename A>
				void record(const CrystalStructure<A>&, const double elapsedTime);

				template <typename A>
				void record(const CrystalStructure<A>&, const size_type elapsedCount, const std::string& evaluationValueName, const double evaluationValue);

				template <typename A>
				void record(const CrystalStructure<A>&, const double elapsedTime, const std::string& evaluationValueName, const double evaluationValue);

				void release() noexcept override;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void writeNumTotalAtom(const size_type numTotalAtom, StreamWriter&) const;

				void writeElapsedTime(const size_type elapsedCount, StreamWriter&) const;
				void writeElapsedTime(const double elapsedTime, StreamWriter&) const;
				void writeEvaluationValue(const std::string& name, const double value, StreamWriter&) const;

				void writeUnitCell(const NumericalMatrix&, StreamWriter&) const;
				void writeUnitCell(const UnitCell&, StreamWriter&) const;

				template <typename A>
				void writeAtomicArrangement(const std::vector<A>&, StreamWriter&) const;

			// Public methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private utility

				void writeCartesianCoordinate(const NumericalVector&, StreamWriter&) const;

			// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				mutable size_type m_frameNumber;
				mutable size_type m_elapsedCount;
				mutable double m_elapsedTime;
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

inline std::string DftToolkit::Openmx::IO::OpenmxMdStreamWriter::correctExtension() noexcept
{
	return std::string{ ".md" };
}

inline bool DftToolkit::Openmx::IO::OpenmxMdStreamWriter::isCorrectExtension(const std::filesystem::path& filePath)
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

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystallographicStructure<A>& structure, const size_type elapsedCount)
{
	if (m_elapsedCount < elapsedCount)
	{
		m_elapsedCount = elapsedCount;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(structure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedCount, streamWriter);
			writeUnitCell(structure.unitCell(), streamWriter);
			writeAtomicArrangement(structure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed count is invalid." };
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystallographicStructure<A>& structure, const double elapsedTime)
{
	if (m_elapsedTime < elapsedTime)
	{
		m_elapsedTime = elapsedTime;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(structure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedTime, streamWriter);
			writeUnitCell(structure.unitCell().basisVectors(), streamWriter);
			writeAtomicArrangement(structure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed time is invalid." };
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystalStructure<A>& crystalStructure, const size_type elapsedCount)
{
	if (m_elapsedCount < elapsedCount)
	{
		m_elapsedCount = elapsedCount;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(crystalStructure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedCount, streamWriter);
			writeUnitCell(crystalStructure.unitCell().basisVectors(), streamWriter);
			writeAtomicArrangement(crystalStructure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed count is invalid." };
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystalStructure<A>& crystalStructure, const double elapsedTime)
{
	if (m_elapsedTime < elapsedTime)
	{
		m_elapsedTime = elapsedTime;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(crystalStructure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedTime, streamWriter);
			writeUnitCell(crystalStructure.basisVectors(), streamWriter);
			writeAtomicArrangement(crystalStructure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed time is invalid." };
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystalStructure<A>& crystalStructure, const size_type elapsedCount, const std::string& evaluationValueName, const double evaluationValue)
{
	if (m_elapsedCount < elapsedCount)
	{
		m_elapsedCount = elapsedCount;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(crystalStructure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedCount, streamWriter);
			writeEvaluationValue(evaluationValueName, evaluationValue, streamWriter);
			writeUnitCell(crystalStructure.basisVectors(), streamWriter);
			writeAtomicArrangement(crystalStructure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed count is invalid." };
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::record(const CrystalStructure<A>& crystalStructure, const double elapsedTime, const std::string& evaluationValueName, const double evaluationValue)
{
	if (m_elapsedTime < elapsedTime)
	{
		m_elapsedTime = elapsedTime;

		StreamWriter streamWriter;
		{
			writeNumTotalAtom(crystalStructure.atoms().size(), streamWriter);
			writeElapsedTime(elapsedTime, streamWriter);
			writeEvaluationValue(evaluationValueName, evaluationValue, streamWriter);
			writeUnitCell(crystalStructure.basisVectors(), streamWriter);
			writeAtomicArrangement(crystalStructure.atoms(), streamWriter);
		}

		FileStream::write(streamWriter.allLines());
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "record", "Argument elapsed time is invalid." };
}

inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::release() noexcept
{
	m_frameNumber = 0;
	m_elapsedCount = 0;
	m_elapsedTime = 0.0;

	FileStream::release();
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

inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::writeNumTotalAtom(const size_type numTotalAtom, StreamWriter& streamWriter) const
{
	streamWriter.writeLine(numTotalAtom);
}

inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::writeElapsedTime(const size_type elapsedCount, StreamWriter& streamWriter) const
{
	streamWriter.write("  elapsed_count=  ");
	streamWriter.write(elapsedCount);
	streamWriter.write(" ()  ");
}

inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::writeUnitCell(const UnitCell& unitCell, StreamWriter& streamWriter) const
{
	writeUnitCell(unitCell.basisVectors(), streamWriter);
}

inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::writeEvaluationValue(const std::string& name, const double value, StreamWriter& streamWriter) const
{
	streamWriter.write(name);
	streamWriter.write("=  ");
	streamWriter.write(value);
	streamWriter.write("  ");
}

template <typename A>
inline void DftToolkit::Openmx::IO::OpenmxMdStreamWriter::writeAtomicArrangement(const std::vector<A>& atomicArrangement, StreamWriter& streamWriter) const
{
	for (const auto& atom : atomicArrangement)
	{
		streamWriter.write("  ");
		streamWriter.write(atom.atomicNumber().toElementSymbol());
		streamWriter.write("  ");

		writeCartesianCoordinate(atom.cartesianCoordinate(), streamWriter);

		streamWriter.write("  ");
		streamWriter.write("0.0  0.0  0.0");
		streamWriter.write("  ");
		streamWriter.write("0.0  0.0  0.0");
		streamWriter.write("  ");
		streamWriter.writeLine("0.0  0.0  0.0  0.0");
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


#endif // !DFTTOOLKIT_OPENMX_IO_OPENMXMDSTREAMWRITER_H
