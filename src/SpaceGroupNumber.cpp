#include "SpaceGroupNumber.h"

#include <limits>

using namespace ChemToolkit::Crystallography::Symmetry;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

SpaceGroupNumber::SpaceGroupNumber() noexcept
	: _number{ 1 }
{
}

SpaceGroupNumber::SpaceGroupNumber(const size_type number)
	: _number{ number }
{
	if (!(isValidNumber(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is out of definitions." };
}

SpaceGroupNumber::SpaceGroupNumber(const std::size_t number)
	: _number{ static_cast<size_type>(number) }
{
	if (std::numeric_limits<size_type>::max() < number)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is too large." };

	if (!(isValidNumber(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is out of definitions." };
}

SpaceGroupNumber::SpaceGroupNumber(const int number)
	: _number{ static_cast<size_type>(number) }
{
	if (std::numeric_limits<size_type>::max() < number)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is too large." };

	if(number < 1)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is less than one." };

	if (!(isValidNumber(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is out of definitions." };
}

SpaceGroupNumber::SpaceGroupNumber(const std::string& symbol)
	: _number{ SpaceGroupSymbol{ symbol }.toSpaceGroupNumber() }
{
}

SpaceGroupNumber::SpaceGroupNumber(const SpaceGroupSymbol& symbol)
	: _number{ symbol.toSpaceGroupNumber() }
{
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

SpaceGroupSymbol SpaceGroupNumber::toSpaceGroupSymbol() const
{
	switch (_number)
	{
	case 1:
		return SpaceGroupSymbol{ std::string{ "P1" } };
	case 2:
		return SpaceGroupSymbol{ std::string{ "P-1" } };
	case 3:
		return SpaceGroupSymbol{ std::string{ "P2" } };
	case 4:
		return SpaceGroupSymbol{ std::string{ "P2_1" } };
	case 5:
		return SpaceGroupSymbol{ std::string{ "C2" } };
	case 6:
		return SpaceGroupSymbol{ std::string{ "Pm" } };
	case 7:
		return SpaceGroupSymbol{ std::string{ "Pc" } };
	case 8:
		return SpaceGroupSymbol{ std::string{ "Cm" } };
	case 9:
		return SpaceGroupSymbol{ std::string{ "Cc" } };
	case 10:
		return SpaceGroupSymbol{ std::string{ "P2/m" } };
	case 11:
		return SpaceGroupSymbol{ std::string{ "P2_1/m" } };
	case 12:
		return SpaceGroupSymbol{ std::string{ "C2/m" } };
	case 13:
		return SpaceGroupSymbol{ std::string{ "P2/c" } };
	case 14:
		return SpaceGroupSymbol{ std::string{ "P2_1/c" } };
	case 15:
		return SpaceGroupSymbol{ std::string{ "C2/c" } };
	case 16:
		return SpaceGroupSymbol{ std::string{ "P222" } };
	case 17:
		return SpaceGroupSymbol{ std::string{ "P222_1" } };
	case 18:
		return SpaceGroupSymbol{ std::string{ "P2_12_12" } };
	case 19:
		return SpaceGroupSymbol{ std::string{ "P2_12_12_1" } };
	case 20:
		return SpaceGroupSymbol{ std::string{ "C222_1" } };
	case 21:
		return SpaceGroupSymbol{ std::string{ "C222" } };
	case 22:
		return SpaceGroupSymbol{ std::string{ "F222" } };
	case 23:
		return SpaceGroupSymbol{ std::string{ "I222" } };
	case 24:
		return SpaceGroupSymbol{ std::string{ "I2_12_12_1" } };
	case 25:
		return SpaceGroupSymbol{ std::string{ "Pmm2" } };
	case 26:
		return SpaceGroupSymbol{ std::string{ "Pmc2_1" } };
	case 27:
		return SpaceGroupSymbol{ std::string{ "Pcc2" } };
	case 28:
		return SpaceGroupSymbol{ std::string{ "Pma2" } };
	case 29:
		return SpaceGroupSymbol{ std::string{ "Pca2_1" } };
	case 30:
		return SpaceGroupSymbol{ std::string{ "Pnc2" } };
	case 31:
		return SpaceGroupSymbol{ std::string{ "Pmn2_1" } };
	case 32:
		return SpaceGroupSymbol{ std::string{ "Pba2" } };
	case 33:
		return SpaceGroupSymbol{ std::string{ "Pna2_1" } };
	case 34:
		return SpaceGroupSymbol{ std::string{ "Pnn2" } };
	case 35:
		return SpaceGroupSymbol{ std::string{ "Cmm2" } };
	case 36:
		return SpaceGroupSymbol{ std::string{ "Cmc2_1" } };
	case 37:
		return SpaceGroupSymbol{ std::string{ "Ccc2" } };
	case 38:
		return SpaceGroupSymbol{ std::string{ "Amm2" } };
	case 39:
		return SpaceGroupSymbol{ std::string{ "Aem2" } };
	case 40:
		return SpaceGroupSymbol{ std::string{ "Ama2" } };
	case 41:
		return SpaceGroupSymbol{ std::string{ "Aea2" } };
	case 42:
		return SpaceGroupSymbol{ std::string{ "Fmm2" } };
	case 43:
		return SpaceGroupSymbol{ std::string{ "Fdd2" } };
	case 44:
		return SpaceGroupSymbol{ std::string{ "Imm2" } };
	case 45:
		return SpaceGroupSymbol{ std::string{ "Iba2" } };
	case 46:
		return SpaceGroupSymbol{ std::string{ "Ima2" } };
	case 47:
		return SpaceGroupSymbol{ std::string{ "Pmmm" } };
	case 48:
		return SpaceGroupSymbol{ std::string{ "Pnnn" } };
	case 49:
		return SpaceGroupSymbol{ std::string{ "Pccm" } };
	case 50:
		return SpaceGroupSymbol{ std::string{ "Pban" } };
	case 51:
		return SpaceGroupSymbol{ std::string{ "Pmma" } };
	case 52:
		return SpaceGroupSymbol{ std::string{ "Pnna" } };
	case 53:
		return SpaceGroupSymbol{ std::string{ "Pmna" } };
	case 54:
		return SpaceGroupSymbol{ std::string{ "Pcca" } };
	case 55:
		return SpaceGroupSymbol{ std::string{ "Pbam" } };
	case 56:
		return SpaceGroupSymbol{ std::string{ "Pccn" } };
	case 57:
		return SpaceGroupSymbol{ std::string{ "Pbcm" } };
	case 58:
		return SpaceGroupSymbol{ std::string{ "Pnnm" } };
	case 59:
		return SpaceGroupSymbol{ std::string{ "Pmmn" } };
	case 60:
		return SpaceGroupSymbol{ std::string{ "Pbcn" } };
	case 61:
		return SpaceGroupSymbol{ std::string{ "Pbca" } };
	case 62:
		return SpaceGroupSymbol{ std::string{ "Pnma" } };
	case 63:
		return SpaceGroupSymbol{ std::string{ "Cmcm" } };
	case 64:
		return SpaceGroupSymbol{ std::string{ "Cmce" } };
	case 65:
		return SpaceGroupSymbol{ std::string{ "Cmmm" } };
	case 66:
		return SpaceGroupSymbol{ std::string{ "Cccm" } };
	case 67:
		return SpaceGroupSymbol{ std::string{ "Cmme" } };
	case 68:
		return SpaceGroupSymbol{ std::string{ "Ccce" } };
	case 69:
		return SpaceGroupSymbol{ std::string{ "Fmmm" } };
	case 70:
		return SpaceGroupSymbol{ std::string{ "Fddd" } };
	case 71:
		return SpaceGroupSymbol{ std::string{ "Immm" } };
	case 72:
		return SpaceGroupSymbol{ std::string{ "Ibam" } };
	case 73:
		return SpaceGroupSymbol{ std::string{ "Ibca" } };
	case 74:
		return SpaceGroupSymbol{ std::string{ "Imma" } };
	case 75:
		return SpaceGroupSymbol{ std::string{ "P4" } };
	case 76:
		return SpaceGroupSymbol{ std::string{ "P4_1" } };
	case 77:
		return SpaceGroupSymbol{ std::string{ "P4_2" } };
	case 78:
		return SpaceGroupSymbol{ std::string{ "P4_3" } };
	case 79:
		return SpaceGroupSymbol{ std::string{ "I4" } };
	case 80:
		return SpaceGroupSymbol{ std::string{ "I4_1" } };
	case 81:
		return SpaceGroupSymbol{ std::string{ "P-4" } };
	case 82:
		return SpaceGroupSymbol{ std::string{ "I-4" } };
	case 83:
		return SpaceGroupSymbol{ std::string{ "P4/m" } };
	case 84:
		return SpaceGroupSymbol{ std::string{ "P4_2/m" } };
	case 85:
		return SpaceGroupSymbol{ std::string{ "P4/n" } };
	case 86:
		return SpaceGroupSymbol{ std::string{ "P4_2/n" } };
	case 87:
		return SpaceGroupSymbol{ std::string{ "I4/m" } };
	case 88:
		return SpaceGroupSymbol{ std::string{ "I4_1/a" } };
	case 89:
		return SpaceGroupSymbol{ std::string{ "P442" } };
	case 90:
		return SpaceGroupSymbol{ std::string{ "P42_12" } };
	case 91:
		return SpaceGroupSymbol{ std::string{ "P4_122" } };
	case 92:
		return SpaceGroupSymbol{ std::string{ "P4_12_12" } };
	case 93:
		return SpaceGroupSymbol{ std::string{ "P4_222" } };
	case 94:
		return SpaceGroupSymbol{ std::string{ "P4_22_12" } };
	case 95:
		return SpaceGroupSymbol{ std::string{ "P4_322" } };
	case 96:
		return SpaceGroupSymbol{ std::string{ "P4_32_12" } };
	case 97:
		return SpaceGroupSymbol{ std::string{ "I422" } };
	case 98:
		return SpaceGroupSymbol{ std::string{ "I4_122" } };
	case 99:
		return SpaceGroupSymbol{ std::string{ "P4mm" } };
	case 100:
		return SpaceGroupSymbol{ std::string{ "P4bm" } };
	case 101:
		return SpaceGroupSymbol{ std::string{ "P4_2cm" } };
	case 102:
		return SpaceGroupSymbol{ std::string{ "P4_2nm" } };
	case 103:
		return SpaceGroupSymbol{ std::string{ "P4cc" } };
	case 104:
		return SpaceGroupSymbol{ std::string{ "P4nc" } };
	case 105:
		return SpaceGroupSymbol{ std::string{ "P4_2mc" } };
	case 106:
		return SpaceGroupSymbol{ std::string{ "P4_2bc" } };
	case 107:
		return SpaceGroupSymbol{ std::string{ "I4mm" } };
	case 108:
		return SpaceGroupSymbol{ std::string{ "I4cm" } };
	case 109:
		return SpaceGroupSymbol{ std::string{ "I4_1md" } };
	case 110:
		return SpaceGroupSymbol{ std::string{ "I4_1cd" } };
	case 111:
		return SpaceGroupSymbol{ std::string{ "P42m" } };
	case 112:
		return SpaceGroupSymbol{ std::string{ "P-42c" } };
	case 113:
		return SpaceGroupSymbol{ std::string{ "P42_1m" } };
	case 114:
		return SpaceGroupSymbol{ std::string{ "P42_1c" } };
	case 115:
		return SpaceGroupSymbol{ std::string{ "P-4m2" } };
	case 116:
		return SpaceGroupSymbol{ std::string{ "P-4c2" } };
	case 117:
		return SpaceGroupSymbol{ std::string{ "P-4b2" } };
	case 118:
		return SpaceGroupSymbol{ std::string{ "P-4n2" } };
	case 119:
		return SpaceGroupSymbol{ std::string{ "I-4m2" } };
	case 120:
		return SpaceGroupSymbol{ std::string{ "I-4c2" } };
	case 121:
		return SpaceGroupSymbol{ std::string{ "I-42m" } };
	case 122:
		return SpaceGroupSymbol{ std::string{ "I-42d" } };
	case 123:
		return SpaceGroupSymbol{ std::string{ "P4/mmm" } };
	case 124:
		return SpaceGroupSymbol{ std::string{ "P4/mcc" } };
	case 125:
		return SpaceGroupSymbol{ std::string{ "P4/nbm" } };
	case 126:
		return SpaceGroupSymbol{ std::string{ "P4/nnc" } };
	case 127:
		return SpaceGroupSymbol{ std::string{ "P4/mbm" } };
	case 128:
		return SpaceGroupSymbol{ std::string{ "P4/mmc" } };
	case 129:
		return SpaceGroupSymbol{ std::string{ "P4/nmm" } };
	case 130:
		return SpaceGroupSymbol{ std::string{ "P4/ncc" } };
	case 131:
		return SpaceGroupSymbol{ std::string{ "P4_2/mmc" } };
	case 132:
		return SpaceGroupSymbol{ std::string{ "P4_2/mcm" } };
	case 133:
		return SpaceGroupSymbol{ std::string{ "P4_2/nbc" } };
	case 134:
		return SpaceGroupSymbol{ std::string{ "P4_2/nnm" } };
	case 135:
		return SpaceGroupSymbol{ std::string{ "P4_2/mbc" } };
	case 136:
		return SpaceGroupSymbol{ std::string{ "P4_2/mnm" } };
	case 137:
		return SpaceGroupSymbol{ std::string{ "P4_2/nmc" } };
	case 138:
		return SpaceGroupSymbol{ std::string{ "P4_2/ncm" } };
	case 139:
		return SpaceGroupSymbol{ std::string{ "I4/mmm" } };
	case 140:
		return SpaceGroupSymbol{ std::string{ "I4/mcm" } };
	case 141:
		return SpaceGroupSymbol{ std::string{ "I4_1/amd" } };
	case 142:
		return SpaceGroupSymbol{ std::string{ "I4_1/acd" } };
	case 143:
		return SpaceGroupSymbol{ std::string{ "P3" } };
	case 144:
		return SpaceGroupSymbol{ std::string{ "P3_1" } };
	case 145:
		return SpaceGroupSymbol{ std::string{ "P3_2" } };
	case 146:
		return SpaceGroupSymbol{ std::string{ "R3" } };
	case 147:
		return SpaceGroupSymbol{ std::string{ "P-3" } };
	case 148:
		return SpaceGroupSymbol{ std::string{ "R-3" } };
	case 149:
		return SpaceGroupSymbol{ std::string{ "P312" } };
	case 150:
		return SpaceGroupSymbol{ std::string{ "P321" } };
	case 151:
		return SpaceGroupSymbol{ std::string{ "P3_112" } };
	case 152:
		return SpaceGroupSymbol{ std::string{ "P3_121" } };
	case 153:
		return SpaceGroupSymbol{ std::string{ "P3_212" } };
	case 154:
		return SpaceGroupSymbol{ std::string{ "P3_221" } };
	case 155:
		return SpaceGroupSymbol{ std::string{ "R32" } };
	case 156:
		return SpaceGroupSymbol{ std::string{ "P3m1" } };
	case 157:
		return SpaceGroupSymbol{ std::string{ "P31m" } };
	case 158:
		return SpaceGroupSymbol{ std::string{ "P3c1" } };
	case 159:
		return SpaceGroupSymbol{ std::string{ "P31c" } };
	case 160:
		return SpaceGroupSymbol{ std::string{ "R3m" } };
	case 161:
		return SpaceGroupSymbol{ std::string{ "R3c" } };
	case 162:
		return SpaceGroupSymbol{ std::string{ "P-31m" } };
	case 163:
		return SpaceGroupSymbol{ std::string{ "P-31c" } };
	case 164:
		return SpaceGroupSymbol{ std::string{ "P-3m1" } };
	case 165:
		return SpaceGroupSymbol{ std::string{ "P-3c1" } };
	case 166:
		return SpaceGroupSymbol{ std::string{ "R-3m" } };
	case 167:
		return SpaceGroupSymbol{ std::string{ "R-3c" } };
	case 168:
		return SpaceGroupSymbol{ std::string{ "P6" } };
	case 169:
		return SpaceGroupSymbol{ std::string{ "P6_1" } };
	case 170:
		return SpaceGroupSymbol{ std::string{ "P6_5" } };
	case 171:
		return SpaceGroupSymbol{ std::string{ "P6_2" } };
	case 172:
		return SpaceGroupSymbol{ std::string{ "P6_4" } };
	case 173:
		return SpaceGroupSymbol{ std::string{ "P6_3" } };
	case 174:
		return SpaceGroupSymbol{ std::string{ "P-6" } };
	case 175:
		return SpaceGroupSymbol{ std::string{ "P6/m" } };
	case 176:
		return SpaceGroupSymbol{ std::string{ "P6_3/m" } };
	case 177:
		return SpaceGroupSymbol{ std::string{ "P622" } };
	case 178:
		return SpaceGroupSymbol{ std::string{ "P6_122" } };
	case 179:
		return SpaceGroupSymbol{ std::string{ "P6_522" } };
	case 180:
		return SpaceGroupSymbol{ std::string{ "P6_222" } };
	case 181:
		return SpaceGroupSymbol{ std::string{ "P6_422" } };
	case 182:
		return SpaceGroupSymbol{ std::string{ "P6_322" } };
	case 183:
		return SpaceGroupSymbol{ std::string{ "P6mm" } };
	case 184:
		return SpaceGroupSymbol{ std::string{ "P6cc" } };
	case 185:
		return SpaceGroupSymbol{ std::string{ "P6_3cm" } };
	case 186:
		return SpaceGroupSymbol{ std::string{ "P6_3mc" } };
	case 187:
		return SpaceGroupSymbol{ std::string{ "P-6m2" } };
	case 188:
		return SpaceGroupSymbol{ std::string{ "P-6c2" } };
	case 189:
		return SpaceGroupSymbol{ std::string{ "P-62m" } };
	case 190:
		return SpaceGroupSymbol{ std::string{ "P-62c" } };
	case 191:
		return SpaceGroupSymbol{ std::string{ "P6/mmm" } };
	case 192:
		return SpaceGroupSymbol{ std::string{ "P6/mcc" } };
	case 193:
		return SpaceGroupSymbol{ std::string{ "P6_3/mcm" } };
	case 194:
		return SpaceGroupSymbol{ std::string{ "P6_3/mmc" } };
	case 195:
		return SpaceGroupSymbol{ std::string{ "P23" } };
	case 196:
		return SpaceGroupSymbol{ std::string{ "F23" } };
	case 197:
		return SpaceGroupSymbol{ std::string{ "I23" } };
	case 198:
		return SpaceGroupSymbol{ std::string{ "P2_13" } };
	case 199:
		return SpaceGroupSymbol{ std::string{ "I2_13" } };
	case 200:
		return SpaceGroupSymbol{ std::string{ "Pm-3" } };
	case 201:
		return SpaceGroupSymbol{ std::string{ "Pn-3" } };
	case 202:
		return SpaceGroupSymbol{ std::string{ "Fm-3" } };
	case 203:
		return SpaceGroupSymbol{ std::string{ "Pd-3" } };
	case 204:
		return SpaceGroupSymbol{ std::string{ "Im-3" } };
	case 205:
		return SpaceGroupSymbol{ std::string{ "Pa-3" } };
	case 206:
		return SpaceGroupSymbol{ std::string{ "Ia-3" } };
	case 207:
		return SpaceGroupSymbol{ std::string{ "P432" } };
	case 208:
		return SpaceGroupSymbol{ std::string{ "P4_232" } };
	case 209:
		return SpaceGroupSymbol{ std::string{ "F432" } };
	case 210:
		return SpaceGroupSymbol{ std::string{ "F4_132" } };
	case 211:
		return SpaceGroupSymbol{ std::string{ "I432" } };
	case 212:
		return SpaceGroupSymbol{ std::string{ "P4_332" } };
	case 213:
		return SpaceGroupSymbol{ std::string{ "P4_132" } };
	case 214:
		return SpaceGroupSymbol{ std::string{ "I4_132" } };
	case 215:
		return SpaceGroupSymbol{ std::string{ "P-43m" } };
	case 216:
		return SpaceGroupSymbol{ std::string{ "F-43m" } };
	case 217:
		return SpaceGroupSymbol{ std::string{ "I-43m" } };
	case 218:
		return SpaceGroupSymbol{ std::string{ "P-43n" } };
	case 219:
		return SpaceGroupSymbol{ std::string{ "F-43c" } };
	case 220:
		return SpaceGroupSymbol{ std::string{ "I-43d" } };
	case 221:
		return SpaceGroupSymbol{ std::string{ "Pm-3m" } };
	case 222:
		return SpaceGroupSymbol{ std::string{ "Pn-3n" } };
	case 223:
		return SpaceGroupSymbol{ std::string{ "Pm-3n" } };
	case 224:
		return SpaceGroupSymbol{ std::string{ "Pn-3m" } };
	case 225:
		return SpaceGroupSymbol{ std::string{ "Fm-3m" } };
	case 226:
		return SpaceGroupSymbol{ std::string{ "Fm-3c" } };
	case 227:
		return SpaceGroupSymbol{ std::string{ "Fd-3m" } };
	case 228:
		return SpaceGroupSymbol{ std::string{ "Fd-3c" } };
	case 229:
		return SpaceGroupSymbol{ std::string{ "Im-3m" } };
	case 230:
		return SpaceGroupSymbol{ std::string{ "Ia-3d" } };
	default:
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "toSpaceGroupSymbol", "Space group number is out of definitions." };
	}
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
