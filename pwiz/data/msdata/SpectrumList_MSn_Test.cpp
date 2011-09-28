//
// $Id$
//
//
// Original author: Barbara Frewen <frewen@u.washington.edu>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#include "SpectrumList_MSn.hpp"
#include "Serializer_MSn.hpp"
#include "TextWriter.hpp"
#include "pwiz/utility/minimxml/XMLWriter.hpp"
#include "pwiz/utility/misc/unit.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Base64.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace pwiz::msdata;
using namespace pwiz::util;
using namespace pwiz::minimxml;


ostream* os_ = 0;

const char *testMSn = 
"H	CreationDate	10-6-2007\n"
"H	Extractor	MakeMS2\n"
"H	ExtractorVersion	2.0\n"
"H	Comments	MakeMS2 written by Michael J. MacCoss, Michael R. Hoopmann, 2007\n"
"H	ExtractorOptions	MS2\n"
"S	116	116	536.39\n"
"Z	2	1071.77\n"
"Z	3	1607.15\n"
"I	RTime	0.4462\n"
"175.4067 0.6\n"
"195.9831 0.9\n"
"236.2524 0.7\n"
"266.1989 0.8\n"
"276.3145 0.6\n"
"278.6099 0.7\n"
"281.1050 62.6\n"
"298.4888 9.2\n"
"299.1033 4\n"
"303.5076 1.8\n"
"330.1379 1.9\n"
"337.1581 0.6\n"
"341.0460 1.8\n"
"342.3756 0.8\n"
"359.0190 1.7\n"
"363.0510 1.1\n"
"373.2335 1.7\n"
"377.2952 1\n"
"399.2092 1.6\n"
"399.8853 1\n"
"403.3747 2.2\n"
"405.3385 0.6\n"
"408.1845 0.8\n"
"409.8174 0.8\n"
"414.0231 1.2\n"
"415.0792 1.8\n"
"417.1754 1.6\n"
"419.3268 1.2\n"
"420.4443 1.5\n"
"421.3249 1.1\n"
"429.1662 2.5\n"
"435.1794 0.8\n"
"436.3020 1.9\n"
"439.8971 1\n"
"440.6216 1.1\n"
"443.7952 0.6\n"
"444.5042 1.2\n"
"447.0575 0.6\n"
"448.5712 1.8\n"
"451.1549 2\n"
"452.1009 0.8\n"
"453.0457 0.7\n"
"461.8017 0.7\n"
"464.1340 2.7\n"
"469.8256 1.8\n"
"471.8412 1\n"
"473.1831 3.2\n"
"474.2579 1\n"
"479.8830 1.1\n"
"482.2438 1.7\n"
"483.1778 2.8\n"
"483.8499 0.8\n"
"486.1272 2.8\n"
"487.5856 0.7\n"
"489.5583 4.5\n"
"490.6985 4.2\n"
"491.4770 5.2\n"
"492.3076 7.9\n"
"496.1183 1.2\n"
"498.5404 2.6\n"
"500.5744 7.8\n"
"501.2284 6.6\n"
"501.9958 3.9\n"
"503.0387 37.4\n"
"505.0599 2.9\n"
"507.2243 1.2\n"
"508.2107 0.8\n"
"509.1219 2.4\n"
"510.4375 0.9\n"
"511.3940 0.7\n"
"513.4534 3.5\n"
"514.4604 4.5\n"
"517.0839 6.3\n"
"518.1203 18.8\n"
"519.0070 65\n"
"525.1636 3.1\n"
"526.4990 9.4\n"
"527.5856 12.2\n"
"560.5077 1.2\n"
"585.4477 2\n"
"621.8749 0.7\n"
"632.6031 1\n"
"634.1111 0.7\n"
"636.3641 1.2\n"
"638.4987 1.1\n"
"640.5447 1.7\n"
"650.5433 0.9\n"
"664.1221 0.8\n"
"698.4615 0.9\n"
"709.7639 1.1\n"
"711.2064 1.1\n"
"726.0311 1.2\n"
"740.0786 0.8\n"
"745.3728 0.6\n"
"757.8849 1\n"
"761.9862 0.8\n"
"774.4131 1.3\n"
"788.2714 1\n"
"829.2268 0.7\n"
"840.0249 2\n"
"856.4430 0.8\n"
"857.6420 0.7\n"
"898.4391 0.8\n"
"902.4149 0.7\n"
"942.4218 1.4\n"
"1026.6023 1\n"
"S	118	118	464.98\n"
"Z	2	928.95\n"
"Z	3	1392.92\n"
"I	RTime	0.4573\n"
"159.1265 1\n"
"176.0755 1.4\n"
"189.2380 1.3\n"
"199.4232 0.8\n"
"205.2997 0.6\n"
"213.3207 1.4\n"
"221.2078 0.9\n"
"231.0154 1.6\n"
"238.9865 1\n"
"244.2399 1.4\n"
"249.2524 1.1\n"
"252.1188 0.9\n"
"253.2228 0.9\n"
"263.1987 0.8\n"
"269.3484 0.8\n"
"272.3980 0.9\n"
"273.2263 0.9\n"
"276.9475 2.2\n"
"279.2698 1.4\n"
"299.1149 0.8\n"
"302.9478 0.8\n"
"305.6744 1.4\n"
"308.3992 1.2\n"
"317.9594 1.3\n"
"322.2583 1.2\n"
"332.9634 0.7\n"
"337.3151 0.8\n"
"346.9624 1.3\n"
"349.1566 2.3\n"
"351.1241 1.3\n"
"357.0767 1.8\n"
"361.8666 0.9\n"
"363.1213 2.1\n"
"365.3057 1.3\n"
"370.7582 2.3\n"
"375.0994 1.5\n"
"377.1262 4.7\n"
"383.2674 2.4\n"
"385.3204 1.1\n"
"386.8699 0.8\n"
"390.2235 2\n"
"391.1315 1.8\n"
"393.2798 1.3\n"
"399.0021 0.7\n"
"399.8439 1\n"
"401.0514 2\n"
"401.9986 2.1\n"
"405.0557 2.3\n"
"406.2401 1.8\n"
"407.1901 2.9\n"
"409.2251 0.9\n"
"410.0710 1.1\n"
"411.0744 4.5\n"
"412.0661 3.5\n"
"416.1904 0.8\n"
"418.3093 0.8\n"
"419.3824 4.1\n"
"420.0984 3.3\n"
"421.1917 10.7\n"
"422.2375 2.3\n"
"423.5126 1.9\n"
"424.3980 1.3\n"
"426.1027 0.9\n"
"426.8067 0.8\n"
"428.1773 5.9\n"
"429.0244 4.4\n"
"430.1233 1.5\n"
"433.7563 1.5\n"
"435.2169 1.1\n"
"435.9414 0.9\n"
"437.1042 2.3\n"
"438.0438 1.2\n"
"442.4830 4.6\n"
"443.4222 1.6\n"
"446.2477 23.2\n"
"447.0900 44.6\n"
"447.9651 5.1\n"
"452.9314 2.1\n"
"455.4204 6.1\n"
"456.3987 3\n"
"464.9924 12\n"
"466.0859 1.3\n"
"482.7336 1.8\n"
"531.3920 0.8\n"
"646.3096 1.8\n"
;

const char *testBMS2 = 
"AwAAAAIAAABDcmVhdGlvbkRhdGUJMTAtNi0yMDA3CgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEV4dHJhY3RvcglNYWtlTVMyCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARXh0cmFjdG9yVmVyc2lvbgkyLjAKAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABDb21tZW50cwlNYWtlTVMyIHdyaXR0ZW4gYnkgTWljaGFlbCBKLiBNYWNDb3NzLCBNaWNoYWVsIFIuIEhvb3BtYW5uLCAyMDA3CgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEV4dHJhY3Rvck9wdGlvbnMJTVMyCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADgjlkAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDOfwAAAAAAAAAAAAAAAAAAAAAAAAAAACGcUwAAAAAA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";


const char *testCMS2 = 
"BAAAAAIAAABDcmVhdGlvbkRhdGUJMTAtNi0yMDA3CgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEV4dHJhY3RvcglNYWtlTVMyCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARXh0cmFjdG9yVmVyc2lvbgkyLjAKAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABDb21tZW50cwlNYWtlTVMyIHdyaXR0ZW4gYnkgTWljaGFlbCBKLiBNYWNDb3NzLCBNaWNoYWVsIFIuIEhvb3BtYW5uLCAyMDA3CgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEV4dHJhY3Rvck9wdGlvbnMJTVMyCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADgjlkAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDOfwAAAAAAAAAAAAAAAAAAAAAAAAAAACGcUwAAAAAA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";

const char* testMS2_v3 = 
"H       CreationDate    11-20-2008\n"
"H       Extractor       MakeMS2\n"
"H       ExtractorVersion        2.0\n"
"H       Comments        MakeMS2 written by Michael J. MacCoss, Michael R. Hoopmann, 2007\n"
"H       ExtractorOptions        MS1/MS2\n"
"S       36      36      612.1900\n"
"I       RTime   0.4338\n"
"I       BPM     0.0000\n"
"I       ConvA   0.0000\n"
"I       TIC     0.00\n"
"I       IIT     0.0000\n"
"I       EZ      1       611.1855        0.8897  27106.6\n"
"Z       1       611.1855\n"
"221.0247 4.1\n"
"222.2416 3\n"
"223.1555 5.9\n"
"240.0833 0\n"
"249.9288 0\n"
"251.0432 2\n"
"252.9677 1.2\n"
"265.0794 0\n"
"267.1087 12.3\n"
"268.0749 6\n"
"269.0090 11.6\n"
"281.1318 7.6\n"
"282.0841 2.1\n"
"283.1004 5.8\n"
"285.1190 1.7\n"
"286.0294 1.8\n"
"287.0404 4.4\n"
"297.2826 0\n"
"298.6739 0\n"
"300.0451 0\n"
"342.9260 0\n"
"347.1423 0\n"
"355.0954 72.4\n"
"356.0631 45\n"
"357.0933 109.3\n"
"S       508     508     441.2300\n"
"I       RTime   6.2752\n"
"I       ConvA   0.0000\n"
"I       ConvB   0.0000\n"
"I       TIC     0.00\n"
"I       IIT     0.0000\n"
"I       EZ      3       1318.7270       5.9467  28240.3\n"
"I       EZ      2       880.4527        6.2403  34674.3\n"
"Z       3       1318.7270\n"
"Z       2       880.4527\n"
"129.1926 60.4\n"
"130.0556 6\n"
"131.1174 1.4\n"
"136.2554 0\n"
"137.9548 2.1\n"
"138.8827 0\n"
"140.1849 1.4\n"
"141.1434 0\n"
"143.1377 2.2\n"
"147.1954 11.1\n"
"150.1817 0\n"
"151.1815 1.5\n"
"152.0610 1.2\n"
"155.0635 1.2\n"
"156.1205 2.6\n"
"157.1824 0\n"
"158.3834 2.3\n"
"161.1998 0\n"
"166.0671 1\n"
"169.2280 1.2\n"
"170.5823 0\n"
"173.1148 0\n"
"175.1335 11.9\n"
"176.1112 1.1\n"
"178.1093 1.3\n"
"180.1308 1.1\n"
"181.2527 1.5\n"
"182.9846 2.5\n"
;

const char* testCMS2_v3 =
"BAAAAAMAAAAgICAgICBDcmVhdGlvbkRhdGUgICAgMTEtMjAtMjAwOAoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAgICAgIEV4dHJhY3RvciAgICAgICBNYWtlTVMyCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAICAgICAgRXh0cmFjdG9yVmVyc2lvbiAgICAgICAgMi4wCgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAgICBDb21tZW50cyAgICAgICAgTWFrZU1TMiB3cml0dGVuIGJ5IE1pY2hhZWwgSi4gTWFjQ29zcywgTWljaGFlbCBSLiBIb29wbWFubiwgMjAwNwoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAgICAgIEV4dHJhY3Rvck9wdGlvbnMgICAgICAgIE1TMS9NUzIKAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD21UoAAAAAACQAAAAkAAAA7FG4HoUhg0AJG94+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAABAAAAGQAAAAEAAABEi2znexmDQAEAAABEi2znexmDQGHDYz8zxdNG0wAAAE4AAAB4nAHIADf/woanV8qga0DWVuwvu8drQOXQItv55GtAOPjCZKoCbkDLEMe6uD1vQP5l9+RhYW9A24r9Zfeeb0DXNO84RZFwQDQRNjy9sXBAh6dXyjLBcECgGi/dJNBwQM4ZUdobknFAf9k9eVihcUBLyAc9m7FxQJZDi2zn0XFACmgibHjgcUAibHh6pfBxQJjdk4eFlHJAZMxdS8iqckDLEMe6uMByQIlBYOXQbnVAiGNd3EaydUCdgCbChjF2QAtGJXUCQXZAeAskKH5RdkCwlF86eJwzNm52YGBwcDh7Zg+QhgOHWTNn2oMYZ8+4ODIwHADyLR2NjT87pKWxAdk7gfimfVraM/uzZ3qQ9cHB2TMTnBgYTJxmzbzlBAB8CRv9/AEAAPwBAABI4XoUrpN7QHDOyEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAAIAAAAcAAAAAwAAAMUgsHLompRAAgAAADJ3LSGfg4tAAwAAAMUgsHLompRAXku+QJqg3EYCAAAAMnctIZ+Di0CKsMdATXIHR+sAAABRAAAAeJwB4AAf/7WmeccpJmBAC7WmecdBYECdEaW9wWNgQMDsnjwsCGFAEce6uI0+YUC8BRIUP1xhQPmgZ7PqhWFA48eYu5akYUAYldQJaORhQG40gLdAZmJA3bWEfNDFYkArhxbZzuViQMuhRbbzAWNAeekmMQhiY0Bg5dAi24NjQEtZhjjWpWNAKqkT0ETMY0C1N/jCZCZkQPkx5q4lwmRAnu+nxksnZUC/DpwzolJlQJayDHGso2VAg8DKoUXkZUAWak3zjgNmQH2utmJ/Q2ZAI0p7gy+EZkBgdk8eFqhmQIj029eB32ZAwPlnynicmzWz0ImB4YCDsfFmewYgSEtjcwDRMP7ZMzwOs2YaOjKAwQH7WTNngnFamhpUnbADRK4BLM4ABWlpdo5nz/QA1S2zB9EgvQwMCg4ANGkdIA==";

void test(SpectrumListPtr sl) 
{
    if (os_)
    {
        TextWriter write(*os_);
        write(*sl);
        *os_ << endl;
    }

    // check easy functions

    unit_assert(sl.get());
    unit_assert(sl->size() == 2);
    unit_assert(sl->find("scan=116") == 0);
    unit_assert(sl->find("scan=118") == 1);

    // check scan 0

    unit_assert(sl->spectrumIdentity(0).index == 0);
    unit_assert(sl->spectrumIdentity(0).id == "scan=116");
    unit_assert(sl->spectrumIdentity(0).sourceFilePosition != -1);

    SpectrumPtr s = sl->spectrum(0, false);

    unit_assert(s.get());
    unit_assert(s->id == "scan=116");
    unit_assert(s->index == 0);
    unit_assert(s->sourceFilePosition != -1);
    unit_assert(s->cvParam(MS_ms_level).valueAs<int>() == 2);
    unit_assert_equal(s->cvParam(MS_total_ion_current).valueAs<double>(), 385.4, 5e-1);
    unit_assert_equal(s->cvParam(MS_base_peak_intensity).valueAs<double>(), 65.0, 5e-1);

    unit_assert(s->precursors.size() == 1);
    Precursor& precursor0 = s->precursors[0];
    unit_assert(precursor0.selectedIons.size() == 1);
    unit_assert_equal(precursor0.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>(), 536.39, 5e-2);

    // This test spectrum only has possible charge states (2 total, values 2,3)
    unit_assert(precursor0.selectedIons[0].cvParam(MS_charge_state).value.empty());
    vector<string> charges;
    BOOST_FOREACH(CVParam& param, precursor0.selectedIons[0].cvParams)
    {
        if (param.cvid == MS_possible_charge_state)
        {
            charges.push_back(param.value);
        }
    }
    unit_assert(charges.size() == 2);
    vector<string>::const_iterator charge_it = charges.begin();
    unit_assert(*charge_it == "2");			
    unit_assert(*(++charge_it) == "3");
    
    unit_assert(s->defaultArrayLength == 106);
    unit_assert(s->binaryDataArrayPtrs.size() == 2);
    unit_assert(s->binaryDataArrayPtrs[0]->hasCVParam(MS_m_z_array));
    unit_assert(s->binaryDataArrayPtrs[1]->hasCVParam(MS_intensity_array));
    unit_assert(s->binaryDataArrayPtrs[0]->data.empty() && s->binaryDataArrayPtrs[1]->data.empty());

    s = sl->spectrum(0, true);
    unit_assert(s->defaultArrayLength == 106);
    unit_assert(s->binaryDataArrayPtrs.size() == 2);
    unit_assert(!s->binaryDataArrayPtrs[0]->data.empty() && !s->binaryDataArrayPtrs[1]->data.empty());

    vector<MZIntensityPair> pairs;
    s->getMZIntensityPairs(pairs);

    if (os_)
    {
        *os_ << "scan 0:\n";
        copy(pairs.begin(), pairs.end(), ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }


    // check scan 1

    unit_assert(sl->spectrumIdentity(1).index == 1);
    unit_assert(sl->spectrumIdentity(1).id == "scan=118");

    s = sl->spectrum(1, true);
    unit_assert(s.get());
    unit_assert(s->id == "scan=118");
    unit_assert(s->index == 1);
    unit_assert(s->sourceFilePosition != -1);
    unit_assert(s->cvParam(MS_ms_level).valueAs<int>() == 2);
    unit_assert(s->scanList.scans.size() == 1);
    unit_assert_equal(s->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds(), 0.4573*60, 5e-4);

    unit_assert(s->precursors.size() == 1);
    Precursor& precursor1 = s->precursors[0];
    unit_assert(precursor1.selectedIons.size() == 1);
    unit_assert_equal(precursor1.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>(), 464.98, 1e-5);

    // This test spectrum only has possible charge states (2 total, values 2,3)
    unit_assert(precursor1.selectedIons[0].cvParam(MS_charge_state).value.empty());
    charges.clear();
    BOOST_FOREACH(CVParam& param, precursor1.selectedIons[0].cvParams)
    {
        if (param.cvid == MS_possible_charge_state)
        {
            charges.push_back(param.value);
        }
    }
    unit_assert(charges.size() == 2);
    charge_it = charges.begin();
    unit_assert(*charge_it == "2");			
    unit_assert(*(++charge_it) == "3");

    unit_assert(s->defaultArrayLength == 85);

    pairs.clear();
    s->getMZIntensityPairs(pairs);

    unit_assert(s->defaultArrayLength == pairs.size());

    if (os_)
    {
        *os_ << "scan 1:\n";
        copy(pairs.begin(), pairs.end(), ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }
}

void test_v3(SpectrumListPtr sl) 
{
    if (os_)
    {
        TextWriter write(*os_);
        write(*sl);
        *os_ << endl;
    }

    // check easy functions

    unit_assert(sl.get());
    unit_assert(sl->size() == 2);
    unit_assert(sl->find("scan=36") == 0);
    unit_assert(sl->find("scan=508") == 1);

    // check scan 0

    unit_assert(sl->spectrumIdentity(0).index == 0);
    unit_assert(sl->spectrumIdentity(0).id == "scan=36");
    unit_assert(sl->spectrumIdentity(0).sourceFilePosition != -1);

    SpectrumPtr s = sl->spectrum(0, false);

    unit_assert(s.get());
    unit_assert(s->id == "scan=36");
    unit_assert(s->index == 0);
    unit_assert(s->sourceFilePosition != -1);
    unit_assert(s->cvParam(MS_ms_level).valueAs<int>() == 2);
    unit_assert_equal(s->cvParam(MS_total_ion_current).valueAs<double>(), 296.2, 5e-1);
    unit_assert_equal(s->cvParam(MS_base_peak_intensity).valueAs<double>(), 109.3, 5e-1);

    unit_assert(s->precursors.size() == 1);
    Precursor& precursor0 = s->precursors[0];
    unit_assert(precursor0.selectedIons.size() == 1);
    unit_assert_equal(precursor0.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>(), 612.19, 5e-2);

    // This test spectrum one charge state, so no possible charge states
    unit_assert(precursor0.selectedIons[0].cvParam(MS_possible_charge_state).value.empty());
    vector<string> charges;
    vector<double> masses;
    BOOST_FOREACH(CVParam& param, precursor0.selectedIons[0].cvParams)
    {
        if (param.cvid == MS_charge_state)
        {
            charges.push_back(param.value);
        }
        if (param.cvid == MS_accurate_mass)
        {
          masses.push_back(lexical_cast<double>(param.value));
        }
    }
    unit_assert(charges.size() == 1);
    vector<string>::const_iterator charge_it = charges.begin();
    unit_assert(*charge_it == "1");			

    unit_assert(masses.size() == 1);
    vector<double>::const_iterator mass_it = masses.begin();
    unit_assert_equal(*mass_it, 611.1855, 5e-4);			
    
    unit_assert(s->defaultArrayLength == 25);
    unit_assert(s->binaryDataArrayPtrs.size() == 2);
    unit_assert(s->binaryDataArrayPtrs[0]->hasCVParam(MS_m_z_array));
    unit_assert(s->binaryDataArrayPtrs[1]->hasCVParam(MS_intensity_array));
    unit_assert(s->binaryDataArrayPtrs[0]->data.empty() && s->binaryDataArrayPtrs[1]->data.empty());

    s = sl->spectrum(0, true);
    unit_assert(s->defaultArrayLength == 25);
    unit_assert(s->binaryDataArrayPtrs.size() == 2);
    unit_assert(!s->binaryDataArrayPtrs[0]->data.empty() && !s->binaryDataArrayPtrs[1]->data.empty());

    vector<MZIntensityPair> pairs;
    s->getMZIntensityPairs(pairs);

    if (os_)
    {
        *os_ << "scan 0:\n";
        copy(pairs.begin(), pairs.end(), ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }


    // check scan 1

    unit_assert(sl->spectrumIdentity(1).index == 1);
    unit_assert(sl->spectrumIdentity(1).id == "scan=508");

    s = sl->spectrum(1, true);
    unit_assert(s.get());
    unit_assert(s->id == "scan=508");
    unit_assert(s->index == 1);
    unit_assert(s->sourceFilePosition != -1);
    unit_assert(s->cvParam(MS_ms_level).valueAs<int>() == 2);
    unit_assert(s->scanList.scans.size() == 1);
    unit_assert_equal(s->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds()/60, 6.2752, 5e-4);

    unit_assert(s->precursors.size() == 1);
    Precursor& precursor1 = s->precursors[0];
    unit_assert(precursor1.selectedIons.size() == 2);
    unit_assert_equal(precursor1.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>(), 441.23, 1e-5);

    // This test spectrum has two charge states, both known so no possible charges
    unit_assert(precursor1.selectedIons[0].cvParam(MS_possible_charge_state).value.empty());
    charges.clear();
    masses.clear();
    BOOST_FOREACH(SelectedIon& si, precursor1.selectedIons)
    {
      BOOST_FOREACH(CVParam& param, si.cvParams)
      {
        if (param.cvid == MS_charge_state)
        {
            charges.push_back(param.value);
        }
        if (param.cvid == MS_accurate_mass)
        {
          masses.push_back(lexical_cast<double>(param.value));
        }
      }
    }
    unit_assert(charges.size() == 2);
    charge_it = charges.begin();
    unit_assert(*charge_it == "3");			
    unit_assert(*(++charge_it) == "2");

    unit_assert(masses.size() == 2);
    mass_it = masses.begin();
    unit_assert_equal(*mass_it, 1318.7270, 5e-4);			
    unit_assert_equal(*(++mass_it), 880.4527, 5e-4);			
    
    unit_assert(s->defaultArrayLength == 28);

    pairs.clear();
    s->getMZIntensityPairs(pairs);

    unit_assert(s->defaultArrayLength == pairs.size());

    if (os_)
    {
        *os_ << "scan 1:\n";
        copy(pairs.begin(), pairs.end(), ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }
}

int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        
        // dummy would normally be read in from file
        MSData dummy;
        dummy.instrumentConfigurationPtrs.push_back(InstrumentConfigurationPtr(new InstrumentConfiguration("LCQDeca")));
        dummy.instrumentConfigurationPtrs.back()->cvParams.push_back(MS_LCQ_Deca);
        dummy.instrumentConfigurationPtrs.back()->userParams.push_back(UserParam("doobie", "420"));
        
        // TEST MS2 TEXT FORMAT
        if (os_)
        {
            *os_ << "test()\n";
            *os_ << "msn:\n" << testMSn << endl;
        }

        shared_ptr<istream> isMS2(new istringstream(testMSn));
        SpectrumListPtr slMS2 = SpectrumList_MSn::create(isMS2, dummy, MSn_Type_MS2);
        test(slMS2);

        // TEST BMS2 BINARY FORMAT 
        const string& testBMS2Base64Str = testBMS2;
        vector<char> binaryBufferBMS2;
        binaryBufferBMS2.resize(Base64::textToBinarySize(testBMS2Base64Str.size()) + 1, '\0');
        Base64::textToBinary(testBMS2Base64Str.c_str(), testBMS2Base64Str.size(), &binaryBufferBMS2[0]);
        string binaryStringBMS2;
        for (size_t i = 0; i < binaryBufferBMS2.size(); i++)
        {
            binaryStringBMS2 += binaryBufferBMS2[i];
        }

        shared_ptr<istream> isBMS2(new istringstream(binaryStringBMS2));
        SpectrumListPtr slBMS2 = SpectrumList_MSn::create(isBMS2, dummy, MSn_Type_BMS2);
        test(slBMS2);

        // TEST CMS2 COMPRESSED BINARY FORMAT (version 2)
        const string& testCMS2Base64Str = testCMS2;
        vector<char> binaryBufferCMS2;
        binaryBufferCMS2.resize(Base64::textToBinarySize(testCMS2Base64Str.size()) + 1, '\0');
        Base64::textToBinary(testCMS2Base64Str.c_str(), testCMS2Base64Str.size(), &binaryBufferCMS2[0]);
        string binaryStringCMS2;
        for (size_t i = 0; i < binaryBufferCMS2.size(); i++)
        {
            binaryStringCMS2 += binaryBufferCMS2[i];
        }

        shared_ptr<istream> isCMS2(new istringstream(binaryStringCMS2));
        SpectrumListPtr slCMS2 = SpectrumList_MSn::create(isCMS2, dummy, MSn_Type_CMS2);
        test(slCMS2);

        // TEST MS2 TEXT FORMAT WITH EZ LINES (version 3)
        shared_ptr<istream> isMS2_v3(new istringstream(testMS2_v3));
        SpectrumListPtr slMS2_v3 = SpectrumList_MSn::create(isMS2_v3, dummy, MSn_Type_MS2);
        test_v3(slMS2_v3);

        // TEST CMS2 COMPRESSED BINARY FORMAT (version 3)
        const string& testCMS2_v3Base64Str = testCMS2_v3;
        vector<char> binaryBufferCMS2_v3;
        binaryBufferCMS2_v3.resize(Base64::textToBinarySize(testCMS2_v3Base64Str.size()) + 1, '\0');
        Base64::textToBinary(testCMS2_v3Base64Str.c_str(), testCMS2_v3Base64Str.size(), &binaryBufferCMS2_v3[0]);
        string binaryStringCMS2_v3;
        for (size_t i = 0; i < binaryBufferCMS2_v3.size(); i++)
        {
            binaryStringCMS2_v3 += binaryBufferCMS2_v3[i];
        }

        shared_ptr<istream> isCMS2_v3(new istringstream(binaryStringCMS2_v3));
        SpectrumListPtr slCMS2_v3 = SpectrumList_MSn::create(isCMS2_v3, dummy, MSn_Type_CMS2);
        test_v3(slCMS2_v3);

        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
    }
}


