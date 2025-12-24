#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALPRODUCTIONREPORTPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALPRODUCTIONREPORTPARAMETERS_H

#include <filesystem>

#include "StreamReader.h"

#include "CrystalDesignRecordParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Diagnostics
		{
			class CrystalProductionReportParameters
			{
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				CrystalProductionReportParameters() noexcept;
				virtual ~CrystalProductionReportParameters() = default;

				CrystalProductionReportParameters(const CrystalProductionReportParameters&) = default;
				CrystalProductionReportParameters(CrystalProductionReportParameters&&) noexcept = default;
				CrystalProductionReportParameters& operator=(const CrystalProductionReportParameters&) = default;
				CrystalProductionReportParameters& operator=(CrystalProductionReportParameters&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const std::filesystem::path& designedCrystalsDirectoryPath() const noexcept;
				bool needPiSymmetryCrystalData() const noexcept;
				bool needInfeasibleCrystalData() const noexcept;
				bool needExceptionalCrystalData() const noexcept;

				double spaceGroupPrecision() const noexcept;
				const CrystalDesignRecordParameters& crystalDesignRecordParameters() const noexcept;

				static double defaultSpaceGroupPrecision() noexcept;


				void setDesignedCrystalsDirectoryPath(const std::filesystem::path&);
				void setPiSymmetryCrystalDataNecessity(const bool) noexcept;
				void setInfeasibleCrystalDataNecessity(const bool) noexcept;
				void setExceptionalCrystalDataNecessity(const bool) noexcept;

				void setSpaceGroupPrecision(const double);
				void setCrystalDesignRecordParameters(const CrystalDesignRecordParameters&) noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void initialize() noexcept;
				void initialize(const System::IO::StreamReader&);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				bool toDataNecessity(const std::string&) const;

				void validateInitializedValues() const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				std::filesystem::path _designedCrystalsDirectoryPath;
				bool _needPiSymmetryCrystalData;
				bool _needInfeasibleCrystalData;
				bool _needExceptionalCrystalData;

				double _spaceGroupPrecision;
				CrystalDesignRecordParameters _crystalDesignRecordParameters;

				static double s_defaultSpaceGroupPrecision;
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

inline const std::filesystem::path& MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::designedCrystalsDirectoryPath() const noexcept
{
	return _designedCrystalsDirectoryPath;
}

inline bool MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::needPiSymmetryCrystalData() const noexcept
{
	return _needPiSymmetryCrystalData;
}

inline bool MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::needInfeasibleCrystalData() const noexcept
{
	return _needInfeasibleCrystalData;
}

inline bool MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::needExceptionalCrystalData() const noexcept
{
	return _needExceptionalCrystalData;
}

inline double MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::spaceGroupPrecision() const noexcept
{
	return _spaceGroupPrecision;
}

inline const MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecordParameters& MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::crystalDesignRecordParameters() const noexcept
{
	return _crystalDesignRecordParameters;
}

inline double MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::defaultSpaceGroupPrecision() noexcept
{
	return s_defaultSpaceGroupPrecision;
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setDesignedCrystalsDirectoryPath(const std::filesystem::path& directoryPath)
{
	if (directoryPath.empty())
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setDesignedCrystalsDirectoryPath", "Argument path is empty." };
	else
		_designedCrystalsDirectoryPath = directoryPath;
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setPiSymmetryCrystalDataNecessity(const bool necessity) noexcept
{
	_needPiSymmetryCrystalData = necessity;
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setInfeasibleCrystalDataNecessity(const bool necessity) noexcept
{
	_needInfeasibleCrystalData = necessity;
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setExceptionalCrystalDataNecessity(const bool necessity) noexcept
{
	_needExceptionalCrystalData = necessity;
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setSpaceGroupPrecision(const double val)
{
	if (0.0 < val)
		_spaceGroupPrecision = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setSpaceGroupPrecision", "Argument value not more than zero." };
}

inline void MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters::setCrystalDesignRecordParameters(const CrystalDesignRecordParameters& parameters) noexcept
{
	_crystalDesignRecordParameters = parameters;
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_DIAGNOSTICS_CRYSTALPRODUCTIONREPORTPARAMETERS_H
