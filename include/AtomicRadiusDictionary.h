#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_ATOMICRADIUSDICTIONARY_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_ATOMICRADIUSDICTIONARY_H

#include <limits>

#include "StreamReader.h"

#include "IonicAtomicNumber.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			class AtomicRadiusDictionary
			{
				using size_type = unsigned short;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;


				struct KeyHasher
				{
					std::size_t operator()(const IonicAtomicNumber&) const noexcept;
				};

				struct KeyEqual
				{
					bool operator()(const IonicAtomicNumber&, const IonicAtomicNumber&) const noexcept;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			private:
				AtomicRadiusDictionary() noexcept = delete;

			public:
				virtual ~AtomicRadiusDictionary() = default;

				AtomicRadiusDictionary(const AtomicRadiusDictionary&) = default;
				AtomicRadiusDictionary(AtomicRadiusDictionary&&) noexcept = default;
				AtomicRadiusDictionary& operator=(const AtomicRadiusDictionary&) = default;
				AtomicRadiusDictionary& operator=(AtomicRadiusDictionary&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				static void initialize() noexcept;
				static void initialize(const System::IO::StreamReader&);
				
				static double getMinimumCovalentRadius(const IonicAtomicNumber&) noexcept;
				static double getMaximumCovalentRadius(const IonicAtomicNumber&) noexcept;
				static double getMinimumIonicRadius(const IonicAtomicNumber&) noexcept;
				static double getMaximumIonicRadius(const IonicAtomicNumber&) noexcept;
				static double getMinimumIonicRepulsionRadius(const IonicAtomicNumber&) noexcept;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				static void validateCovalentRadii(const double minimum, const double maximum);
				static void validateIonicRadii(const double mininum, const double maximum);
				static void validateIonicRepulsionRadii(const double mininum);

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				static std::unordered_map<IonicAtomicNumber, double, KeyHasher, KeyEqual> s_minCovalentRadiusDictionary;
				static std::unordered_map<IonicAtomicNumber, double, KeyHasher, KeyEqual> s_maxCovalentRadiusDictionary;
				static std::unordered_map<IonicAtomicNumber, double, KeyHasher, KeyEqual> s_minIonicRadiusDictionary;
				static std::unordered_map<IonicAtomicNumber, double, KeyHasher, KeyEqual> s_maxIonicRadiusDictionary;
				static std::unordered_map<IonicAtomicNumber, double, KeyHasher, KeyEqual> s_minIonicRepulsionRadiusDictionary;
			};



			inline std::size_t AtomicRadiusDictionary::KeyHasher::operator()(const IonicAtomicNumber& ian) const noexcept
			{
				return ian.getHashCode();
			}

			inline bool AtomicRadiusDictionary::KeyEqual::operator()(const IonicAtomicNumber& iana, const IonicAtomicNumber& ianb) const noexcept
			{
				return (iana == ianb);
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

inline double MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::getMinimumCovalentRadius(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_minCovalentRadiusDictionary.find(ionicAtomicNumber);

	if (iter == s_minCovalentRadiusDictionary.end())
		return 0.0;
	else
		return iter->second;
}

inline double MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::getMaximumCovalentRadius(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_maxCovalentRadiusDictionary.find(ionicAtomicNumber);

	if (iter == s_maxCovalentRadiusDictionary.end())
		return std::numeric_limits<double>::max();
	else
		return iter->second;
}

inline double MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::getMinimumIonicRadius(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_minIonicRadiusDictionary.find(ionicAtomicNumber);

	if (iter == s_minIonicRadiusDictionary.end())
		return 0.0;
	else
		return iter->second;
}

inline double MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::getMaximumIonicRadius(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_maxIonicRadiusDictionary.find(ionicAtomicNumber);

	if (iter == s_maxIonicRadiusDictionary.end())
		return std::numeric_limits<double>::max();
	else
		return iter->second;
}

inline double MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::getMinimumIonicRepulsionRadius(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_minIonicRepulsionRadiusDictionary.find(ionicAtomicNumber);

	if (iter == s_minIonicRepulsionRadiusDictionary.end())
		return 0.0;
	else
		return iter->second;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_ATOMICRADIUSDICTIONARY_H
