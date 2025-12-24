#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALDESIGNRECORDER_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALDESIGNRECORDER_H

#include <string>
#include <mutex>

#include "OpenmxMdStreamWriter.h"

#include "CrystalDesignRecordParameters.h"

#include "ConstrainingCrystalStructure.h"
#include "ObjectiveCrystalStructure.h"

#include "ConstrainingMolecularStructure.h"
#include "ObjectiveMolecularStructure.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Diagnostics
		{
			class CrystalDesignRecorder
			{
				using size_type = std::size_t;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;

				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;

				using ConstrainingCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure;
				using ObjectiveCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure;

				using ConstrainingMolecularStructure = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingMolecularStructure;
				using ObjectiveMolecularStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure;


				struct DynamicAtom
				{
					DynamicAtom(const AtomicNumber, const NumericalVector&) noexcept;
					virtual ~DynamicAtom() = default;

					DynamicAtom(const DynamicAtom&) = default;
					DynamicAtom(DynamicAtom&&) noexcept = default;
					DynamicAtom& operator=(const DynamicAtom&) = default;
					DynamicAtom& operator=(DynamicAtom&&) noexcept = default;

					AtomicNumber atomicNumber() const noexcept;
					const NumericalVector& cartesianCoordinate() const noexcept;
					NumericalVector& cartesianCoordinate() noexcept;


				private:
					AtomicNumber _atomicNumber;
					NumericalVector _cartesianCoordinate;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				CrystalDesignRecorder() noexcept;
				explicit CrystalDesignRecorder(const CrystalDesignRecordParameters&);

				virtual ~CrystalDesignRecorder() = default;

				CrystalDesignRecorder(const CrystalDesignRecorder&) = default;
				CrystalDesignRecorder(CrystalDesignRecorder&&) noexcept = default;
				CrystalDesignRecorder& operator=(const CrystalDesignRecorder&) = default;
				CrystalDesignRecorder& operator=(CrystalDesignRecorder&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const std::filesystem::path& mpiCrystalProductionDirectoryPath() const noexcept;
				bool needMdRecord() const noexcept;

				void setMpiCrystalProductionDirectoryPath(const std::filesystem::path& mpiCrystalProductionDirectoryPath);
				void setCrystalDesignRecordParameters(const CrystalDesignRecordParameters&);
				void setProductionName(const std::string&);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void record(const ObjectiveCrystalStructure&);
				void record(const ObjectiveMolecularStructure&);

				void forceRecord(const ObjectiveCrystalStructure&);
				void forceRecord(const ObjectiveMolecularStructure&);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				CrystalDesignRecordParameters _crystalDesignRecordParameters;

				std::filesystem::path _mpiCrystalProductionDirectoryPath;
				std::string _productionName;
				DftToolkit::Openmx::IO::OpenmxMdStreamWriter _mdStreamWriter;


				mutable size_type m_elapsedCount;
				mutable size_type m_mdRecordIntervalCount;

				static std::mutex s_recordMutex;
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

inline MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::AtomicNumber MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::DynamicAtom::atomicNumber() const noexcept
{
	return _atomicNumber;
}

inline const MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::NumericalVector& MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::DynamicAtom::cartesianCoordinate() const noexcept
{
	return _cartesianCoordinate;
}

inline MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::NumericalVector& MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::DynamicAtom::cartesianCoordinate() noexcept
{
	return _cartesianCoordinate;
}

inline const std::filesystem::path& MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::mpiCrystalProductionDirectoryPath() const noexcept
{
	return _mpiCrystalProductionDirectoryPath;
}

inline bool MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::needMdRecord() const noexcept
{
	return _crystalDesignRecordParameters.needMdRecord();
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder::setCrystalDesignRecordParameters(const CrystalDesignRecordParameters& parameters)
{
	_crystalDesignRecordParameters = parameters;

	_mpiCrystalProductionDirectoryPath.clear();
	_productionName.clear();
	_mdStreamWriter.release();
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALDESIGNRECORDER_H
