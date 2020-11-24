#!/usr/bin/env python3


import linecache  # Reads specific lines
from time import strftime  # Time and date printout
import glob  # Directory listings for file counting
import os.path  # Verify that files or directories exist
import argparse  # Parse commands at execution
import shutil  # Copy files (XDS.INP)

# Parser definitions

parser = argparse.ArgumentParser()

parser.add_argument(
    "-th", help="cell parameters for thaumatin activated", action="store_true"
)
parser.add_argument(
    "-pk", help="cell parameters for proteinase K activated", action="store_true"
)
parser.add_argument(
    "-index",
    help="fraction of indexed reflections required is lowered",
    action="store_true",
)

args = parser.parse_args()

# Version of I23_XDS_INP_gen
I23_v = "r0.5"

# Variable parameters

kappa0_safe_start_list = []
kappa0_safe_end_list = []

minimum_frames = 10  # minimum number of frames for dataset detection
det_z_offset = (
    -0.8
)  # AW 09/03/17                                     # Z axis offset for detector, from distance in cbf header

bckg_deg = 30  # Standard range in degrees for background calculation
bckg_fail_deg = 90  # Range in degrees for background calculation in case kappa shadow cannot be avoided

kappa0_safe_start_list.append(-145)  # Start of kappa safe region 1, in deg
kappa0_safe_end_list.append(-35)  # End of kappa safe region 1, in deg

kappa0_safe_start_list.append(215)  # Start of kappa safe region 2, in deg
kappa0_safe_end_list.append(270)  # End of kappa safe region 2, in deg

index_fraction_def = 0.5
index_fraction_low = 0.3

cbf_l = "cbf_header.log"  # Log file to create for cbf header

# Parser overview

if args.index:
    index_fraction = index_fraction_low
else:
    index_fraction = index_fraction_def


if args.th and args.pk:
    print("\n    ! Cannot have both th and pk cell active. Aborting ! \n")
    quit()

# Functions for finding lines, line numbers and cbf files
def line_find(
    lookup, f_path
):  # Define line_find for identifying and printing lines based on search keywords and given file (f_path)
    with open(f_path, "r", encoding="ISO-8859-1") as MyFile:
        for num, line in enumerate(MyFile, 1):
            if lookup in line:
                line_find_out = linecache.getline(f_path, num).rstrip("\n")
    return line_find_out
    # Returns line which matched input. CANNOT handle multiple matches --> will only print one


def line_n_find(
    lookup, f_path
):  # Define line_find for identifying and printing lines based on search keywords and given file (f_path)
    with open(f_path, "r", encoding="ISO-8859-1") as MyFile:
        for num, line in enumerate(MyFile, 1):
            if lookup in line:
                line_n_find_out = num + 1
    return line_n_find_out


def cbfCounter(
    path, prefix
):  # Searches a given path for number of cbf files with given prefix
    file_number = len(glob.glob1(path, prefix + "*.cbf"))
    return file_number
    # Return number of cbf files


def frame_calculator(oscillation_f, start_angle, number_of_frames, start, end):

    frames_out = []

    frame_start = int((start - start_angle) / oscillation_f)
    frame_end = int((end - start_angle) / oscillation_f)

    if (
        frame_start < 0
        or frame_end > number_of_frames
        or frame_start > number_of_frames
        or frame_end <= 0
    ):
        print("Error - Frame range calculation out of bonds. Aborting")
        quit()
    else:
        frames_out.append(frame_start)
        frames_out.append(frame_end)
    return frames_out


def frame_calculator(oscillation_f, start_angle, number_of_frames, start, end):

    frames_out = []

    frame_start = int((start - start_angle) / oscillation_f)
    frame_end = int((end - start_angle) / oscillation_f)

    if (
        frame_start < 0
        or frame_end > number_of_frames
        or frame_start > number_of_frames
        or frame_end <= 0
    ):
        print("Error - Frame range calculation out of bonds. Aborting")
        quit()
    else:
        frames_out.append(frame_start)
        frames_out.append(frame_end)
    return frames_out


def background_selection(
    start_angle,
    oscillation_f,
    number_of_frames,
    kappa0_safe_start_list,
    kappa0_safe_end_list,
    bckg_deg,
    bckg_fail_deg,
):

    end_angle = start_angle + (oscillation_f * number_of_frames)
    bckg_solution = False

    for index, value in enumerate(kappa0_safe_start_list):

        kappa0_safe_start = kappa0_safe_start_list[index]
        kappa0_safe_end = kappa0_safe_end_list[index]

        if start_angle <= kappa0_safe_start and end_angle >= (
            kappa0_safe_start + bckg_deg
        ):
            bckg_start_angle = kappa0_safe_start
            bckg_end_angle = kappa0_safe_start + bckg_deg

            frames = frame_calculator(
                oscillation_f,
                start_angle,
                number_of_frames,
                bckg_start_angle,
                bckg_end_angle,
            )
            bckg_start = frames[0]
            bckg_end = frames[1]

            bckg_solution = True
            status_message = (
                "Background selection: Using "
                + str(bckg_deg)
                + " deg of data for background calculation ("
                + str(bckg_start_angle)
                + " deg to "
                + str(bckg_end_angle)
                + " deg) in shadow-free zone (assumed to be "
                + str(kappa0_safe_start)
                + " deg to "
                + str(kappa0_safe_end)
                + " deg for kappa = 0 deg)"
            )

        elif (
            start_angle >= kappa0_safe_start
            and start_angle <= (end_angle - bckg_deg)
            and end_angle >= (start_angle + bckg_deg)
            and (start_angle + bckg_deg) <= kappa0_safe_end
        ):
            bckg_start_angle = start_angle
            bckg_end_angle = start_angle + bckg_deg

            frames = frame_calculator(
                oscillation_f,
                start_angle,
                number_of_frames,
                bckg_start_angle,
                bckg_end_angle,
            )
            bckg_start = frames[0]
            bckg_end = frames[1]

            bckg_solution = True
            status_message = (
                "Background selection: Using "
                + str(bckg_deg)
                + " deg of data for background calculation ("
                + str(bckg_start_angle)
                + " deg to "
                + str(bckg_end_angle)
                + " deg) in shadow-free zone (assumed to be "
                + str(kappa0_safe_start)
                + " deg to "
                + str(kappa0_safe_end)
                + " deg for kappa = 0 deg)"
            )

        if bckg_solution == True:
            break

    if bckg_solution == False:
        if (end_angle - start_angle) >= bckg_fail_deg:
            bckg_start_angle = start_angle
            bckg_end_angle = start_angle + bckg_fail_deg

            frames = frame_calculator(
                oscillation_f,
                start_angle,
                number_of_frames,
                bckg_start_angle,
                bckg_end_angle,
            )
            bckg_start = frames[0]
            bckg_end = frames[1]

            status_message = (
                "Background selection: Lack of shadow-free omega range. Using "
                + str(bckg_fail_deg)
                + " deg of data for background ("
                + str(bckg_start_angle)
                + " deg to "
                + str(bckg_end_angle)
                + " deg)"
            )

        else:
            bckg_start = 1
            bckg_end = number_of_frames

            status_message = (
                "Background selection: Fallback condition of all frames for background calculation ("
                + str(start_angle)
                + " deg to "
                + str(end_angle)
                + " deg, total of "
                + str(end_angle - start_angle)
                + " deg)"
            )

    bckg_start = bckg_start + 1

    bckg_frames_output = []

    bckg_frames_output.append(bckg_start)
    bckg_frames_output.append(bckg_end)
    bckg_frames_output.append(status_message)

    return bckg_frames_output


# Startup screen with arguments printout

print(("""
    I23 XDS.INP generator version """ + str(
    I23_v
) + """

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Assumes kappa goniometer and a detector z offset of """ + str(
    det_z_offset
) + """ mm !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""))

if args.th:
    print("    -th    [Thaumatin cell active]")
if args.pk:
    print("    -pk    [Proteinase K cell active]")
if args.index:
    print(("    -index    [Fraction of indexed reflections required is lowered from " + str(
        int(index_fraction_def * 100)
    ) + "% to " + str(
        int(index_fraction_low * 100)
    ) + "%]"))

print("")

# Check for XDS.INP in directory

XDSpINP_exist = os.path.isfile("XDS.INP")

while XDSpINP_exist == True:
    try:
        XDSpINP_ow = (
            input("    [?] Warning XDS.INP already present. Overwrite (y/n)?    ")
        ).lower()
        if XDSpINP_ow == ("y"):
            print("\n    ! Will overwrite !    (Ctrl + C to abort) \n")
            XDSpINP_exist = False
        if XDSpINP_ow == ("n"):
            print("\n    ! Aborting ! \n")
            quit()
    except KeyboardInterrupt:
        print("\n\n    ! Keyboard abort ! \n")
        quit()

# Path input and image detection

found_dataset = False

while found_dataset == False:
    myPath = ""
    while (
        os.path.isdir(myPath) == False
    ):  # Check that input (.cbf) files actually exist, otherwise force new input
        try:
            myPath = input("    Enter data directory: ")
        except KeyboardInterrupt:
            print("\n\n    ! Keyboard abort !\n")
            quit()
        if (
            os.path.isdir(myPath) == False
        ):  # If input file cannot be found, print warning
            print (
                "\n    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n    !! Could not find given path !! \n    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n"
            )

    if myPath[-1] != "/":
        myPath = myPath + "/"

    cbfset = glob.glob(myPath + "/*_00001.cbf")
    dataset_numbers = len(cbfset)
    data_sets = {}
    data_index = 0
    data_sel = 0

    for n in range(0, dataset_numbers):
        prefix_n = cbfset[n].split("/", 100)[-1][:-9]
        frame_numb_n = int(cbfCounter((myPath + "/"), prefix_n))
        if frame_numb_n >= minimum_frames:
            if data_index == 0:
                print(("\n    Available datasets (minimum " + str(
                    minimum_frames
                ) + " frames):"))
                data_index += 1
                print(("    [" + str(
                    data_index
                ) + "]    " + prefix_n + "    (number of frames: " + str(
                    frame_numb_n
                ) + ")"))
                data_sets[data_index] = prefix_n
                n += 1
            else:
                data_index += 1
                print(("    [" + str(
                    data_index
                ) + "]    " + prefix_n + "    (number of frames: " + str(
                    frame_numb_n
                ) + ")"))
                data_sets[data_index] = prefix_n
                n += 1
        else:
            n = n + 1

    if data_index == 0:
        print(("\n    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n    !! No datasets with >= " + str(
            minimum_frames
        ) + " frames found !!\n    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"))

    if data_index != 0:
        found_dataset = True
        print("")

loop_run = 1

while loop_run == 1:
    try:
        data_sel = int(
            input("    Enter dataset index number for XDS.INP generation: ")
        )
    except ValueError:
        data_sel = int("0")
    except KeyboardInterrupt:
        print("\n\n    ! Keyboard abort ! \n")
        quit()
    if data_sel > 0 and data_sel <= data_index:
        loop_run = 0

# Frame name for header readout and data set name
frame_in = myPath + data_sets[data_sel]  # Combination of path and data prefix
frame = frame_in + "00001.cbf"  # First frame in data set
frame_template = (
    myPath + data_sets[data_sel] + "?????.cbf"
)  # Generalised description of data set cbf frames

# Verify that detector is Pilatus 12M

try:
    P12_line = line_n_find("# Detector: PILATUS 12M, S/N 120-0100", frame)
except UnboundLocalError:
    P12_line = 0
if P12_line == 0:
    print("\n    ! Detector not identified as I23 Pilatus 12M. Verify input. Aborting !\n")
    quit()

# Identifying line numbers in cbf file
start_line = (
    line_n_find("_array_data.header_contents", frame)
) + 1  # Find first line with header information
end_line = (
    line_n_find("_array_data.data", frame)
) - 3  # Find last line with header ifnromation

# Time of file generation

time_of_generation = strftime("%Y-%m-%d %H:%M:%S")

# Create log file with header information
try:
    with open(
        cbf_l, "w"
    ) as file:  # Creates log file with location of cbf frame used. Overwrites file if present
        file.write(
            "Frame location:             " + frame + "\n"
        )  # Writes frame location to log file
        file.write(
            "Date and time for reading header:     " + time_of_generation + "\n" + "\n"
        )  # Prints date and time when reading cbf header
        for i in range(start_line, end_line):  # Appends header from cbf file
            file.write(linecache.getline(frame, i))
except Exception:
    print("\n    ! Something went wrong when attempting to write the cbf header log file. Aborting !\n")
    quit()

# Fetch parameters from header log file for XDS.INP
wavelength = float(
    line_find("Wavelength", cbf_l).rstrip("A \n").lstrip("# Wavelength ")
)  # Wavelength value (A), as float
oscillation_f = float(
    line_find("Angle_increment", cbf_l).rstrip("deg. \n").lstrip("# Angle_increment ")
)  # Oscillation per frame (degree), as float
start_angle = float(
    line_find(lookup="Start_angle", f_path=cbf_l)
    .rstrip("deg. \n")
    .lstrip("# Start_angle ")
)  # Start angle, single axis (degree), as float
detector_distance = (
    float(
        line_find("Detector_distance", cbf_l)
        .rstrip("m. \n")
        .lstrip("# Detector_distance ")
    )
    * 1000
    + det_z_offset
)  # Detector distance converted to mm, as float
overload_cutoff = int(
    line_find("Count_cutoff", cbf_l).rstrip("counts \n").lstrip("# Count_cutoff ")
)  # Saturated pixels
number_of_frames = int(
    cbfCounter(myPath, data_sets[data_sel])
)  # Number of cbf frames in data set, as integer

bckg_frames = background_selection(
    start_angle,
    oscillation_f,
    number_of_frames,
    kappa0_safe_start_list,
    kappa0_safe_end_list,
    bckg_deg,
    bckg_fail_deg,
)

try:
    detector_yoff = (
        float(
            line_find("Detector_Voffset", cbf_l)
            .rstrip("m. \n")
            .lstrip("# Detector_Voffset ")
        )
        * 1000
    )  # Detector vertical offset converted to mm, as float
except UnboundLocalError:
    print("\n    ! Warning: Could not find detector y offset in CBF-header !")
    detector_yoff = 0

detector_yoff_check = False
if detector_yoff != 0:
    print("\n    ! Detector vertical offset detected. Beam center defaults are expected to be incorrect !")
    detector_yoff_check = True
while detector_yoff_check == True:
    dety_wrn = (input("    [?] Continue (y/n)?:    ")).lower()
    if dety_wrn == ("y"):
        detector_yoff_check == False
        break
    if dety_wrn == ("n"):
        print("\n    ! Aborting ! \n")
        quit()

# Cell parameters

if args.th:
    cell_printout = "  SPACE_GROUP_NUMBER= 92 \n  UNIT_CELL_CONSTANTS= 58 58 150 90.000 90.000 90.000"

elif args.pk:
    cell_printout = "  SPACE_GROUP_NUMBER= 96 \n  UNIT_CELL_CONSTANTS= 68 68 102 90.000 90.000 90.000"

else:
    cell_printout = ""

try:
    with open("XDS.INP", "w") as f:  # Generate XDS.INP. Overwrites any previous copy
        f.write(
            """!========================================================================================
!=
!= I23_XDS.INP_gen version: """
            + I23_v
            + """
!=
!= This file follows mostly the order of keywords given in the XDS documentation
!= 
!= Date and time of XDS.INP generation: """
            + time_of_generation
            + """
!= 
!= """
            + str(bckg_frames[2])
            + """
!=
!= For more details please check the online documentation at
!= http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_parameters.html
!=
!========================================================================================
!=
!=
!========================================================================================
 
!===================================================
!=========== Job control
  JOB= ALL
! JOB= DEFPIX INTEGRATE CORRECT
! JOB= IDXREF DEFPIX INTEGRATE CORRECT
! JOB= COLSPOT IDXREF
! JOB= CORRECT

! SECONDS=0                               ! Maximum number of seconds to wait until data image must appear
! MAXIMUM_NUMBER_OF_JOBS=4                ! Speeds-up COLSPOT & INTEGRATE on a Linux-cluster
! MAXIMUM_NUMBER_OF_PROCESSORS=8          ! <32;ignored by single cpu version of xds
! TEST=1                                  ! Test flag. 1,2 additional diagnostics and images

!===================================================
!=========== Detector hardware
  DETECTOR= PILATUS
  NX= 2463
  NY= 5071
  QX= 0.172000
  QY= 0.172000
  OVERLOAD= """
            + str(overload_cutoff)
            + """
  MINIMUM_VALID_PIXEL_VALUE= 0
  TRUSTED_REGION= 0.0 1.45
  SENSOR_THICKNESS= 0.320
  VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 1000 40000        ! Verify by inspection of ABS.cbf
 
!===================================================
!=========== Trusted detector region
  INCLUDE_RESOLUTION_RANGE= 999.0 0.0
 
!===================================================
!=========== Detector geometry
  DIRECTION_OF_DETECTOR_X-AXIS= -1  0  0
  DIRECTION_OF_DETECTOR_Y-AXIS=  0 -1  0
  DETECTOR_DISTANCE= """
            + str(detector_distance)
            + """


  ORGX= -153.7                          ! difference between beam on detector and 1230
  ORGY= -9.8


!  ORGX= -152.30                           ! old settings (before 04/04/2019)
!  ORGY= -40.00
!  ORGX= -157.00                           ! old settings (before 01/12/2018)
!  ORGY= -0.57
 
!===================================================
!=========== Data images
  NAME_TEMPLATE_OF_DATA_FRAMES="""
            + frame_template
            + """ CBF
  DATA_RANGE= 1 """
            + str(number_of_frames)
            + """
  BACKGROUND_RANGE= """
            + str(bckg_frames[0])
            + """ """
            + str(bckg_frames[1])
            + """
  SPOT_RANGE= 1 """
            + str(number_of_frames)
            + """
 
!===================================================
!=========== Rotation axis
  ROTATION_AXIS= 1.000000  0.000000  0.000000             ! change "1.000000" to "-1.000000" if day 1 goniometer
  OSCILLATION_RANGE= """
            + str(oscillation_f)
            + """
 
!===================================================
!=========== Incident beam
  X-RAY_WAVELENGTH= """
            + str(wavelength)
            + """
  INCIDENT_BEAM_DIRECTION= -0.001653 -0.000977 0.403273
  FRACTION_OF_POLARIZATION= 0.9999
  POLARIZATION_PLANE_NORMAL= 0 1 0
  AIR= 0.0
 
!===================================================
!=========== Crystal

! SPACE_GROUP_NUMBER= 96                    ! Pk
! UNIT_CELL_CONSTANTS= 68 68 102 90.000 90.000 90.000        ! Pk
! SPACE_GROUP_NUMBER= 92                    ! Th
! UNIT_CELL_CONSTANTS= 58 58 150 90.000 90.000 90.000        ! Th
! SPACE_GROUP_NUMBER= 230                   ! Germanate
! UNIT_CELL_CONSTANTS= 51.3 51.3 51.3 90.0 90.0 90.0               ! Germanate
! SPACE_GROUP_NUMBER= n                    
! UNIT_CELL_CONSTANTS= a b c x y z

"""
            + cell_printout
            + """

  FRIEDEL'S_LAW= FALSE
  STARTING_ANGLE= """
            + str(start_angle)
            + """
  STARTING_FRAME= 1

! REFERENCE_DATA_SET=../path_to/XDS_ASCII.HKL  ! change path and re-run CORRECT
 
!===================================================
!=========== Background and peak pixels
! STRONG_PIXEL= 4                        
  MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 2
 
!===================================================
!========== Beamstop
! UNTRUSTED_RECTANGLE= 3 1030 2570 2610
! UNTRUSTED_ELLIPSE= 1024 1118 2547 2642
!
!========== border issue upper left side
! UNTRUSTED_RECTANGLE= 0 486 423 428
! UNTRUSTED_RECTANGLE= 0 486 615 620
! UNTRUSTED_RECTANGLE= 480 486 424 618
! UNTRUSTED_RECTANGLE= 0 5 424 618
!
!========== noisy chips (25-06-2018)
! UNTRUSTED_RECTANGLE= 0 62 1060 1158
! UNTRUSTED_RECTANGLE= 0 62 1908 2006
!
!========== noisy modules (18-09-2019)
! UNTRUSTED_RECTANGLE= 988 1476 3180 3376
! UNTRUSTED_RECTANGLE= 0 488 3604 3800

!========== noisy modules (26-11-2019)
! UNTRUSTED_RECTANGLE= 0 488 3604 3800

!========== blank chip Y-3 X+2 modules from bs (21-02-20)
! UNTRUSTED_RECTANGLE= 2035 2098 3180 3279

!========== blank module Y+11 X+2 modules from bs (21-02-20)
! UNTRUSTED_RECTANGLE= 1975 2462 211 407

!========== noisy module Y-5 X+2 modules from bs (21-02-20)
! UNTRUSTED_RECTANGLE= 1975 2403 3603 3702

!========== dead module NOV 2020
 UNTRUSTED_RECTANGLE= 1975 2403 3390 3588

!===================================================
!=========== Indexing
  MAXIMUM_ERROR_OF_SPOT_POSITION= 2.0
  MINIMUM_FRACTION_OF_INDEXED_SPOTS= """
            + str(index_fraction)
            + """

!Additional parameters for fine tuning that rarely need to be changed
! INDEX_ERROR=0.05 INDEX_MAGNITUDE=8 INDEX_QUALITY=0.8
! SEPMIN=4.0       ! default is 6 for other detectors
! CLUSTER_RADIUS=2 ! default is 3 for other detectors
! MAXIMUM_ERROR_OF_SPINDLE_POSITION=2.0

!===================================================
!============== Decision constants for finding crystal symmetry
! Constants for detection of lattice symmetry (IDXREF, CORRECT)
! MAX_CELL_AXIS_ERROR=0.03                     ! Maximum relative error in cell axes tolerated 
! MAX_CELL_ANGLE_ERROR=2.0                     ! Maximum cell angle error tolerated

! Constants for detection of space group symmetry (CORRECT) 
! TEST_RESOLUTION_RANGE=8.0 4.5                 ! Resolution range should cover a sufficient number of strong reflections.
! MIN_RFL_Rmeas= 50                         ! Minimum #reflections needed for calculation of Rmeas
! MAX_FAC_Rmeas= 0.0                        ! Is this correct?
 
!===================================================
!=========== Refinement
  REFINE(IDXREF)= POSITION BEAM AXIS ORIENTATION CELL     ! SEGMENT    !/ BEAM AXIS ORIENTATION CELL    (old settings)        / Use all frames for spot finding for these defaults
  REFINE(INTEGRATE)= POSITION BEAM ORIENTATION CELL     ! AXIS        !/ BEAM ORIENTATION CELL     (old settings)        / Can exclude POSITION and check if cell parameters are stable
  REFINE(CORRECT)= BEAM ORIENTATION CELL AXIS POSITION    ! SEGMENT    !/ BEAM AXIS ORIENTATION CELL    (old settings)        / If "low resolution", remove POSITION from CORRECT, INTEGRATE, IDXREF

! DEFAULT_REFINE_SEGMENT= POSITION ORIENTATION
! MINIMUM_NUMBER_OF_REFLECTIONS/SEGMENT=300
 
!===================================================
!===========  Ice rings and exclusions
! EXCLUDE_RESOLUTION_RANGE=   3.911    3.868 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   3.676    3.637 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   3.455    3.415 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   2.686    2.644 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   2.263    2.230 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   2.081    2.051 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   1.960    1.931 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   1.927    1.901 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   1.893    1.867 ! pre-defined ice-ring resolution
! EXCLUDE_RESOLUTION_RANGE=   1.726    1.708 ! pre-defined ice-ring resolution

!Exclude fat-ring
! EXCLUDE_RESOLUTION_RANGE=   6.85 6.65

! MINIMUM_ZETA=0.05                     ! Defines width of 'blind region' (XPLAN,INTEGRATE,CORRECT)
! WFAC1=1.0                          ! Controls number of rejected MISFITS in CORRECT, a larger value leads to fewer rejections.
! REJECT_ALIEN=20.0                     ! Automatic rejection of very strong reflections

!===================================================
!=========== Peak profiles
!Specification of the peak profile parameters below overrides the automatic determination from the images
!Suggested values are listed near the end of INTEGRATE.LP
! BEAM_DIVERGENCE=   0.80                     ! arctan(spot diameter/DETECTOR_DISTANCE)
! BEAM_DIVERGENCE_E.S.D.=   0.080             ! half-width (Sigma) of BEAM_DIVERGENCE
! REFLECTING_RANGE=  0.780                 ! for crossing the Ewald sphere on shortest route
! REFLECTING_RANGE_E.S.D.=  0.113             ! half-width (mosaicity) of REFLECTING_RANGE

  NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA= 15    ! used by: INTEGRATE  / P12M default example "13"
  NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 15         ! used by: INTEGRATE  / P12M default example "9"

! DELPHI= 6.0                        ! controls the number of reference profiles and scaling factors
! CUT= 2.0                            ! defines the integration region for profile fitting
! MINPK= 75.0                         ! minimum required percentage of observed reflection intensity
 
!===================================================
!======= PARAMETERS CONTROLLING CORRECTION FACTORS (used by CORRECT)
! MINIMUM_I/SIGMA= 3.0                     ! minimum intensity/sigma required for scaling reflections
! NBATCH= -1                          ! controls the number of correction factors along image numbers
! REFLECTIONS/CORRECTION_FACTOR= 50               ! minimum #reflections/correction needed
! PATCH_SHUTTER_PROBLEM= TRUE                     ! FALSE is default
! STRICT_ABSORPTION_CORRECTION= TRUE              ! FALSE is default
! CORRECTIONS= DECAY MODULATION ABSORPTION

!===================================================
!=== Detector segment definitions
! refined against 20141016/germ8979

! SEGMENT 1
!
SEGMENT=  1 487 1 195
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00106 0.00008
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00023 -0.14270 0.98977
SEGMENT_ORGX=  1229.49
SEGMENT_ORGY=  97.59
SEGMENT_DISTANCE=  250.22

!
! SEGMENT 2
!
SEGMENT=  495 981 1 195
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00094 0.00014
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00027 -0.14289 0.98974
SEGMENT_ORGX=  1229.78
SEGMENT_ORGY=  97.27
SEGMENT_DISTANCE=  250.23

!
! SEGMENT 3
!
SEGMENT=  989 1475 1 195
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00011 0.00035
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00033 -0.14157 0.98993
SEGMENT_ORGX=  1231.21
SEGMENT_ORGY=  99.28
SEGMENT_DISTANCE=  250.14

!
! SEGMENT 4
!
SEGMENT=  1483 1969 1 195
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00023 -0.00027
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00030 -0.14316 0.98970
SEGMENT_ORGX=  1231.50
SEGMENT_ORGY=  96.76
SEGMENT_DISTANCE=  250.20

!
! SEGMENT 5
!
SEGMENT=  1977 2463 1 195
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00017 -0.00024
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00026 -0.14305 0.98972
SEGMENT_ORGX=  1231.68
SEGMENT_ORGY=  97.00
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 6
!
SEGMENT=  1 487 213 407
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00051 0.00011
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00011 0.00238 1.00000
SEGMENT_ORGX=  1229.91
SEGMENT_ORGY=  309.52
SEGMENT_DISTANCE=  250.07

!
! SEGMENT 7
!
SEGMENT=  495 981 213 407
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00090 -0.00018
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00018 0.00238 1.00000
SEGMENT_ORGX=  1229.56
SEGMENT_ORGY=  309.61
SEGMENT_DISTANCE=  250.17

!
! SEGMENT 8
!
SEGMENT=  989 1475 213 407
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00023 -0.00046
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00046 0.00229 1.00000
SEGMENT_ORGX=  1231.29
SEGMENT_ORGY=  309.33
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 9
!
SEGMENT=  1483 1969 213 407
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00020 -0.00009
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00009 0.00194 1.00000
SEGMENT_ORGX=  1231.24
SEGMENT_ORGY=  308.88
SEGMENT_DISTANCE=  250.16

!
! SEGMENT 10
!
SEGMENT=  1977 2463 213 407
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00028 -0.00048
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00048 0.00221 1.00000
SEGMENT_ORGX=  1230.89
SEGMENT_ORGY=  308.99
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 11
!
SEGMENT=  1 487 425 619
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00021 -0.00034
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00036 0.14700 0.98914
SEGMENT_ORGX=  1230.66
SEGMENT_ORGY=  521.26
SEGMENT_DISTANCE=  249.99

!
! SEGMENT 12
!
SEGMENT=  495 981 425 619
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00098 0.00012
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00002 0.14705 0.98913
SEGMENT_ORGX=  1229.96
SEGMENT_ORGY=  520.71
SEGMENT_DISTANCE=  250.19

!
! SEGMENT 13
!
SEGMENT=  989 1475 425 619
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00034 -0.00060
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00065 0.14764 0.98904
SEGMENT_ORGX=  1230.78
SEGMENT_ORGY=  521.60
SEGMENT_DISTANCE=  250.08

!
! SEGMENT 14
!
SEGMENT=  1483 1969 425 619
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00016 -0.00074
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00076 0.14716 0.98911
SEGMENT_ORGX=  1231.12
SEGMENT_ORGY=  520.91
SEGMENT_DISTANCE=  250.12

!
! SEGMENT 15
!
SEGMENT=  1977 2463 425 619
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00023 -0.00055
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00051 0.14723 0.98910
SEGMENT_ORGX=  1231.96
SEGMENT_ORGY=  520.83
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 16
!
SEGMENT=  1 487 637 831
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00007 0.00022
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00019 0.28893 0.95735
SEGMENT_ORGX=  1230.68
SEGMENT_ORGY=  732.43
SEGMENT_DISTANCE=  250.01

!
! SEGMENT 17
!
SEGMENT=  495 981 637 831
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00070 -0.00020
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00039 0.28888 0.95736
SEGMENT_ORGX=  1230.30
SEGMENT_ORGY=  732.83
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 18
!
SEGMENT=  989 1475 637 831
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00084 -0.00061
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00083 0.28961 0.95715
SEGMENT_ORGX=  1230.42
SEGMENT_ORGY=  733.71
SEGMENT_DISTANCE=  250.09

!
! SEGMENT 19
!
SEGMENT=  1483 1969 637 831
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00004 -0.00059
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00055 0.28951 0.95718
SEGMENT_ORGX=  1231.74
SEGMENT_ORGY=  733.47
SEGMENT_DISTANCE=  250.19

!
! SEGMENT 20
!
SEGMENT=  1977 2463 637 831
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00004 -0.00059
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00058 0.28889 0.95736
SEGMENT_ORGX=  1231.49
SEGMENT_ORGY=  732.31
SEGMENT_DISTANCE=  250.16

!
! SEGMENT 21
!
SEGMENT=  1 487 849 1043
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00014 0.00006
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00001 0.42489 0.90524
SEGMENT_ORGX=  1230.92
SEGMENT_ORGY=  944.90
SEGMENT_DISTANCE=  250.04

!
! SEGMENT 22
!
SEGMENT=  495 981 849 1043
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00076 0.00043
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00007 0.42509 0.90515
SEGMENT_ORGX=  1230.46
SEGMENT_ORGY=  945.23
SEGMENT_DISTANCE=  250.21

!
! SEGMENT 23
!
SEGMENT=  989 1475 849 1043
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00042 -0.00107
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00115 0.42569 0.90487
SEGMENT_ORGX=  1231.58
SEGMENT_ORGY=  946.10
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 24
!
SEGMENT=  1483 1969 849 1043
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00011 -0.00050
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00041 0.42530 0.90505
SEGMENT_ORGX=  1231.79
SEGMENT_ORGY=  945.30
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 25
!
SEGMENT=  1977 2463 849 1043
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00061 -0.00050
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00019 0.42520 0.90510
SEGMENT_ORGX=  1232.66
SEGMENT_ORGY=  945.33
SEGMENT_DISTANCE=  250.38

!
! SEGMENT 26
!
SEGMENT=  1 487 1061 1255
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00024 -0.00062
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00065 0.55176 0.83400
SEGMENT_ORGX=  1231.29
SEGMENT_ORGY=  1157.37
SEGMENT_DISTANCE=  249.97

!
! SEGMENT 27
!
SEGMENT=  495 981 1061 1255
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00104 0.00037
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00026 0.55182 0.83396
SEGMENT_ORGX=  1229.76
SEGMENT_ORGY=  1156.82
SEGMENT_DISTANCE=  250.24

!
! SEGMENT 28
!
SEGMENT=  989 1475 1061 1255
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00016 -0.00047
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00048 0.55207 0.83380
SEGMENT_ORGX=  1231.55
SEGMENT_ORGY=  1157.39
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 29
!
SEGMENT=  1483 1969 1061 1255
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00003 -0.00099
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00084 0.55221 0.83371
SEGMENT_ORGX=  1232.36
SEGMENT_ORGY=  1157.57
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 30
!
SEGMENT=  1977 2463 1061 1255
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00022 -0.00053
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00032 0.55182 0.83396
SEGMENT_ORGX=  1231.97
SEGMENT_ORGY=  1156.86
SEGMENT_DISTANCE=  250.28

!
! SEGMENT 31
!
SEGMENT=  1 487 1273 1467
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00015 0.00013
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00001 0.66732 0.74477
SEGMENT_ORGX=  1231.16
SEGMENT_ORGY=  1369.38
SEGMENT_DISTANCE=  250.07

!
! SEGMENT 32
!
SEGMENT=  495 981 1273 1467
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00063 0.00076
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00015 0.66760 0.74452
SEGMENT_ORGX=  1230.08
SEGMENT_ORGY=  1369.84
SEGMENT_DISTANCE=  250.20

!
! SEGMENT 33
!
SEGMENT=  989 1475 1273 1467
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00011 -0.00026
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00027 0.66758 0.74454
SEGMENT_ORGX=  1231.67
SEGMENT_ORGY=  1370.10
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 34
!
SEGMENT=  1483 1969 1273 1467
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00046 -0.00045
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00003 0.66738 0.74472
SEGMENT_ORGX=  1232.27
SEGMENT_ORGY=  1369.95
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 35
!
SEGMENT=  1977 2463 1273 1467
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00007 -0.00091
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00063 0.66809 0.74408
SEGMENT_ORGX=  1232.06
SEGMENT_ORGY=  1370.62
SEGMENT_DISTANCE=  250.32

!
! SEGMENT 36
!
SEGMENT=  1 487 1485 1679
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00016 0.00111
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00083 0.76825 0.64015
SEGMENT_ORGX=  1230.41
SEGMENT_ORGY=  1580.14
SEGMENT_DISTANCE=  250.22

!
! SEGMENT 37
!
SEGMENT=  495 981 1485 1679
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00094 0.00104
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00005 0.76800 0.64045
SEGMENT_ORGX=  1229.59
SEGMENT_ORGY=  1580.69
SEGMENT_DISTANCE=  250.25

!
! SEGMENT 38
!
SEGMENT=  989 1475 1485 1679
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00049 -0.00114
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00111 0.76871 0.63959
SEGMENT_ORGX=  1232.32
SEGMENT_ORGY=  1582.25
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 39
!
SEGMENT=  1483 1969 1485 1679
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00030 -0.00083
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00076 0.76828 0.64012
SEGMENT_ORGX=  1232.00
SEGMENT_ORGY=  1580.86
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 40
!
SEGMENT=  1977 2463 1485 1679
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00033 -0.00084
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00028 0.76806 0.64038
SEGMENT_ORGX=  1232.42
SEGMENT_ORGY=  1580.84
SEGMENT_DISTANCE=  250.35

!
! SEGMENT 41
!
SEGMENT=  1 487 1697 1891
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00012 0.00028
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00004 0.85314 0.52169
SEGMENT_ORGX=  1230.97
SEGMENT_ORGY=  1793.75
SEGMENT_DISTANCE=  250.14

!
! SEGMENT 42
!
SEGMENT=  495 981 1697 1891
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00052 0.00053
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.85309 0.52176
SEGMENT_ORGX=  1230.50
SEGMENT_ORGY=  1793.76
SEGMENT_DISTANCE=  250.17

!
! SEGMENT 43
!
SEGMENT=  989 1475 1697 1891
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00025 -0.00126
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00045 0.85311 0.52173
SEGMENT_ORGX=  1233.42
SEGMENT_ORGY=  1793.55
SEGMENT_DISTANCE=  250.09

!
! SEGMENT 44
!
SEGMENT=  1483 1969 1697 1891
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00092 -0.00086
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00033 0.85325 0.52150
SEGMENT_ORGX=  1233.35
SEGMENT_ORGY=  1794.29
SEGMENT_DISTANCE=  250.22

!
! SEGMENT 45
!
SEGMENT=  1977 2463 1697 1891
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00008 -0.00113
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00052 0.85326 0.52149
SEGMENT_ORGX=  1232.75
SEGMENT_ORGY=  1793.99
SEGMENT_DISTANCE=  250.40

!
! SEGMENT 46
!
SEGMENT=  1 487 1909 2103
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00028 -0.00004
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00028 0.92003 0.39185
SEGMENT_ORGX=  1231.26
SEGMENT_ORGY=  2006.73
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 47
!
SEGMENT=  495 981 1909 2103
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00049 0.00109
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00002 0.91985 0.39227
SEGMENT_ORGX=  1229.91
SEGMENT_ORGY=  2005.90
SEGMENT_DISTANCE=  250.28

!
! SEGMENT 48
!
SEGMENT=  989 1475 1909 2103
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00073 -0.00102
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00107 0.92009 0.39172
SEGMENT_ORGX=  1232.34
SEGMENT_ORGY=  2007.10
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 49
!
SEGMENT=  1483 1969 1909 2103
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00050 -0.00059
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00023 0.91986 0.39226
SEGMENT_ORGX=  1232.26
SEGMENT_ORGY=  2006.43
SEGMENT_DISTANCE=  250.25

!
! SEGMENT 50
!
SEGMENT=  1977 2463 1909 2103
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00023 -0.00070
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00049 0.92007 0.39176
SEGMENT_ORGX=  1232.60
SEGMENT_ORGY=  2006.89
SEGMENT_DISTANCE=  250.17

!
! SEGMENT 51
!
SEGMENT=  1 487 2121 2315
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00035 0.00033
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00025 0.96695 0.25495
SEGMENT_ORGX=  1230.91
SEGMENT_ORGY=  2217.57
SEGMENT_DISTANCE=  250.21

!
! SEGMENT 52
!
SEGMENT=  495 981 2121 2315
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00057 0.00085
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00033 0.96701 0.25473
SEGMENT_ORGX=  1230.18
SEGMENT_ORGY=  2217.68
SEGMENT_DISTANCE=  250.26

!
! SEGMENT 53
!
SEGMENT=  989 1475 2121 2315
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00024 -0.00056
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00038 0.96692 0.25510
SEGMENT_ORGX=  1232.29
SEGMENT_ORGY=  2217.15
SEGMENT_DISTANCE=  250.12

!
! SEGMENT 54
!
SEGMENT=  1483 1969 2121 2315
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00038 -0.00100
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00062 0.96709 0.25444
SEGMENT_ORGX=  1232.83
SEGMENT_ORGY=  2218.09
SEGMENT_DISTANCE=  250.23

!
! SEGMENT 55
!
SEGMENT=  1977 2463 2121 2315
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00036 -0.00067
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00052 0.96712 0.25433
SEGMENT_ORGX=  1232.46
SEGMENT_ORGY=  2218.04
SEGMENT_DISTANCE=  250.22

!
! SEGMENT 56
!
SEGMENT=  1 487 2333 2527
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00051 -0.00020
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00053 0.99369 0.11217
SEGMENT_ORGX=  1232.01
SEGMENT_ORGY=  2429.32
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 57
!
SEGMENT=  495 981 2333 2527
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00033 0.00099
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00021 0.99374 0.11172
SEGMENT_ORGX=  1230.34
SEGMENT_ORGY=  2429.69
SEGMENT_DISTANCE=  250.24

!
! SEGMENT 58
!
SEGMENT=  989 1475 2333 2527
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00067 -0.00123
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00080 0.99372 0.11190
SEGMENT_ORGX=  1233.56
SEGMENT_ORGY=  2429.39
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 59
!
SEGMENT=  1483 1969 2333 2527
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00007 -0.00038
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00011 0.99374 0.11170
SEGMENT_ORGX=  1232.41
SEGMENT_ORGY=  2429.91
SEGMENT_DISTANCE=  250.15

!
! SEGMENT 60
!
SEGMENT=  1977 2463 2333 2527
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00080 -0.00014
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00081 0.99380 0.11117
SEGMENT_ORGX=  1232.16
SEGMENT_ORGY=  2429.85
SEGMENT_DISTANCE=  250.02

!
! SEGMENT 61
!
SEGMENT=  1 487 2545 2739
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00000 0.00000
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00000 0.99944 -0.03352
SEGMENT_ORGX=  1231.50
SEGMENT_ORGY=  2641.50
SEGMENT_DISTANCE=  250.00

!
! SEGMENT 62
!
SEGMENT=  495 981 2545 2739
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00015 0.00057
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.99942 -0.03399
SEGMENT_ORGX=  1230.79
SEGMENT_ORGY=  2642.03
SEGMENT_DISTANCE=  250.23

!
! SEGMENT 63
!
SEGMENT=  989 1475 2545 2739
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00005 0.00096
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00001 0.99941 -0.03430
SEGMENT_ORGX=  1230.24
SEGMENT_ORGY=  2642.62
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 64
!
SEGMENT=  1483 1969 2545 2739
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00000 0.00000
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00000 0.99944 -0.03352
SEGMENT_ORGX=  1231.50
SEGMENT_ORGY=  2641.50
SEGMENT_DISTANCE=  250.00

!
! SEGMENT 65
!
SEGMENT=  1977 2463 2545 2739
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00000 0.00000
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00000 0.99944 -0.03352
SEGMENT_ORGX=  1231.50
SEGMENT_ORGY=  2641.50
SEGMENT_DISTANCE=  250.00

!
! SEGMENT 66
!
SEGMENT=  1 487 2757 2951
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00049 0.00114
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00068 0.98402 -0.17805
SEGMENT_ORGX=  1230.13
SEGMENT_ORGY=  2853.83
SEGMENT_DISTANCE=  250.28

!
! SEGMENT 67
!
SEGMENT=  495 981 2757 2951
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00038 0.00066
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00049 0.98400 -0.17816
SEGMENT_ORGX=  1230.97
SEGMENT_ORGY=  2853.57
SEGMENT_DISTANCE=  250.19

!
! SEGMENT 68
!
SEGMENT=  989 1475 2757 2951
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00032 -0.00068
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00019 0.98391 -0.17866
SEGMENT_ORGX=  1232.83
SEGMENT_ORGY=  2854.17
SEGMENT_DISTANCE=  250.09

!
! SEGMENT 69
!
SEGMENT=  1483 1969 2757 2951
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00037 -0.00113
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00016 0.98397 -0.17834
SEGMENT_ORGX=  1233.94
SEGMENT_ORGY=  2853.59
SEGMENT_DISTANCE=  250.27

!
! SEGMENT 70
!
SEGMENT=  1977 2463 2757 2951
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00016 -0.00027
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00011 0.98399 -0.17824
SEGMENT_ORGX=  1232.61
SEGMENT_ORGY=  2853.59
SEGMENT_DISTANCE=  250.07

!
! SEGMENT 71
!
SEGMENT=  1 487 2969 3163
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00030 -0.00024
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00021 0.94770 -0.31916
SEGMENT_ORGX=  1232.27
SEGMENT_ORGY=  3065.75
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 72
!
SEGMENT=  495 981 2969 3163
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00012 0.00040
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00024 0.94756 -0.31957
SEGMENT_ORGX=  1231.22
SEGMENT_ORGY=  3066.48
SEGMENT_DISTANCE=  250.19

!
! SEGMENT 73
!
SEGMENT=  989 1475 2969 3163
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00058 -0.00069
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00033 0.94764 -0.31935
SEGMENT_ORGX=  1233.03
SEGMENT_ORGY=  3065.98
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 74
!
SEGMENT=  1483 1969 2969 3163
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00048 -0.00088
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.94769 -0.31918
SEGMENT_ORGX=  1233.23
SEGMENT_ORGY=  3065.67
SEGMENT_DISTANCE=  250.21

!
! SEGMENT 75
!
SEGMENT=  1977 2463 2969 3163
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00165 -0.00010
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00153 0.94753 -0.31967
SEGMENT_ORGX=  1232.01
SEGMENT_ORGY=  3064.97
SEGMENT_DISTANCE=  250.27

!
! SEGMENT 76
!
SEGMENT=  1 487 3181 3375
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00041 -0.00026
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00025 0.89152 -0.45299
SEGMENT_ORGX=  1231.75
SEGMENT_ORGY=  3276.93
SEGMENT_DISTANCE=  250.02

!
! SEGMENT 77
!
SEGMENT=  495 981 3181 3375
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00004 0.00066
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00026 0.89140 -0.45322
SEGMENT_ORGX=  1230.98
SEGMENT_ORGY=  3277.45
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 78
!
SEGMENT=  989 1475 3181 3375
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00007 -0.00097
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00038 0.89142 -0.45318
SEGMENT_ORGX=  1233.06
SEGMENT_ORGY=  3277.26
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 79
!
SEGMENT=  1483 1969 3181 3375
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00042 -0.00120
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00017 0.89131 -0.45339
SEGMENT_ORGX=  1233.60
SEGMENT_ORGY=  3277.98
SEGMENT_DISTANCE=  250.25

!
! SEGMENT 80
!
SEGMENT=  1977 2463 3181 3375
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00012 0.00005
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00013 0.89154 -0.45293
SEGMENT_ORGX=  1231.33
SEGMENT_ORGY=  3276.45
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 81
!
SEGMENT=  1 487 3393 3587
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00041 0.00034
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00053 0.81635 -0.57756
SEGMENT_ORGX=  1231.08
SEGMENT_ORGY=  3489.23
SEGMENT_DISTANCE=  250.08

!
! SEGMENT 82
!
SEGMENT=  495 981 3393 3587
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00031 -0.00018
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00015 0.81619 -0.57779
SEGMENT_ORGX=  1232.06
SEGMENT_ORGY=  3489.30
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 83
!
SEGMENT=  989 1475 3393 3587
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00062 -0.00059
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.81642 -0.57746
SEGMENT_ORGX=  1232.83
SEGMENT_ORGY=  3488.72
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 84
!
SEGMENT=  1483 1969 3393 3587
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00070 -0.00071
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00016 0.81619 -0.57778
SEGMENT_ORGX=  1233.04
SEGMENT_ORGY=  3489.28
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 85
!
SEGMENT=  1977 2463 3393 3587
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00098 -0.00041
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00056 0.81610 -0.57791
SEGMENT_ORGX=  1232.61
SEGMENT_ORGY=  3489.13
SEGMENT_DISTANCE=  250.26

!
! SEGMENT 86
!
SEGMENT=  1 487 3605 3799
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00117 0.00085
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00026 0.72365 -0.69017
SEGMENT_ORGX=  1229.48
SEGMENT_ORGY=  3701.02
SEGMENT_DISTANCE=  250.37

!
! SEGMENT 87
!
SEGMENT=  495 981 3605 3799
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00032 0.00058
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.72347 -0.69035
SEGMENT_ORGX=  1230.16
SEGMENT_ORGY=  3701.80
SEGMENT_DISTANCE=  250.19

!
! SEGMENT 88
!
SEGMENT=  989 1475 3605 3799
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00037 -0.00043
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00003 0.72375 -0.69006
SEGMENT_ORGX=  1232.16
SEGMENT_ORGY=  3701.06
SEGMENT_DISTANCE=  250.08

!
! SEGMENT 89
!
SEGMENT=  1483 1969 3605 3799
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00077 -0.00057
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00017 0.72366 -0.69016
SEGMENT_ORGX=  1233.24
SEGMENT_ORGY=  3701.41
SEGMENT_DISTANCE=  250.18

!
! SEGMENT 90
!
SEGMENT=  1977 2463 3605 3799
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00064 -0.00079
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00008 0.72355 -0.69027
SEGMENT_ORGX=  1232.92
SEGMENT_ORGY=  3701.58
SEGMENT_DISTANCE=  250.26

!
! SEGMENT 91
!
SEGMENT=  1 487 3817 4011
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00020 0.00050
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00052 0.61513 -0.78843
SEGMENT_ORGX=  1230.97
SEGMENT_ORGY=  3915.48
SEGMENT_DISTANCE=  250.16

!
! SEGMENT 92
!
SEGMENT=  495 981 3817 4011
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00045 0.00073
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00030 0.61597 -0.78777
SEGMENT_ORGX=  1230.32
SEGMENT_ORGY=  3913.67
SEGMENT_DISTANCE=  250.20

!
! SEGMENT 93
!
SEGMENT=  989 1475 3817 4011
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00007 -0.00044
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00030 0.61591 -0.78782
SEGMENT_ORGX=  1231.73
SEGMENT_ORGY=  3913.61
SEGMENT_DISTANCE=  250.11

!
! SEGMENT 94
!
SEGMENT=  1483 1969 3817 4011
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00096 0.00050
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00099 0.61532 -0.78828
SEGMENT_ORGX=  1232.34
SEGMENT_ORGY=  3914.15
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 95
!
SEGMENT=  1977 2463 3817 4011
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00022 -0.00030
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00010 0.61556 -0.78809
SEGMENT_ORGX=  1232.00
SEGMENT_ORGY=  3914.38
SEGMENT_DISTANCE=  250.14

!
! SEGMENT 96
!
SEGMENT=  1 487 4029 4223
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00053 0.00024
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00047 0.49448 -0.86919
SEGMENT_ORGX=  1231.71
SEGMENT_ORGY=  4126.96
SEGMENT_DISTANCE=  250.01

!
! SEGMENT 97
!
SEGMENT=  495 981 4029 4223
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00022 0.00073
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00074 0.49423 -0.86933
SEGMENT_ORGX=  1231.30
SEGMENT_ORGY=  4127.33
SEGMENT_DISTANCE=  250.14

!
! SEGMENT 98
!
SEGMENT=  989 1475 4029 4223
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00092 -0.00030
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00020 0.49510 -0.86884
SEGMENT_ORGX=  1233.08
SEGMENT_ORGY=  4125.50
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 99
!
SEGMENT=  1483 1969 4029 4223
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00018 0.00027
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00032 0.49433 -0.86928
SEGMENT_ORGX=  1231.65
SEGMENT_ORGY=  4126.42
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 100
!
SEGMENT=  1977 2463 4029 4223
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00091 -0.00047
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00004 0.49460 -0.86912
SEGMENT_ORGX=  1233.32
SEGMENT_ORGY=  4126.12
SEGMENT_DISTANCE=  250.23

!
! SEGMENT 101
!
SEGMENT=  1 487 4241 4435
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00032 0.00022
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00032 0.36345 -0.93161
SEGMENT_ORGX=  1231.91
SEGMENT_ORGY=  4338.18
SEGMENT_DISTANCE=  250.10

!
! SEGMENT 102
!
SEGMENT=  495 981 4241 4435
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00027 0.00044
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00031 0.36324 -0.93169
SEGMENT_ORGX=  1231.07
SEGMENT_ORGY=  4338.53
SEGMENT_DISTANCE=  250.12

!
! SEGMENT 103
!
SEGMENT=  989 1475 4241 4435
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00101 0.00002
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00038 0.36357 -0.93157
SEGMENT_ORGX=  1232.95
SEGMENT_ORGY=  4337.57
SEGMENT_DISTANCE=  250.05

!
! SEGMENT 104
!
SEGMENT=  1483 1969 4241 4435
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00059 -0.00019
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00003 0.36334 -0.93166
SEGMENT_ORGX=  1232.44
SEGMENT_ORGY=  4337.83
SEGMENT_DISTANCE=  250.09

!
! SEGMENT 105
!
SEGMENT=  1977 2463 4241 4435
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00173 -0.00028
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00036 0.36352 -0.93159
SEGMENT_ORGX=  1234.43
SEGMENT_ORGY=  4337.41
SEGMENT_DISTANCE=  250.26

!
! SEGMENT 106
!
SEGMENT=  1 487 4453 4647
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00058 0.00007
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00020 0.22382 -0.97463
SEGMENT_ORGX=  1232.24
SEGMENT_ORGY=  4551.01
SEGMENT_DISTANCE=  249.95

!
! SEGMENT 107
!
SEGMENT=  495 981 4453 4647
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00071 -0.00041
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00055 0.22448 -0.97448
SEGMENT_ORGX=  1231.01
SEGMENT_ORGY=  4549.95
SEGMENT_DISTANCE=  250.14

!
! SEGMENT 108
!
SEGMENT=  989 1475 4453 4647
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00065 0.00001
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00015 0.22437 -0.97450
SEGMENT_ORGX=  1232.65
SEGMENT_ORGY=  4549.97
SEGMENT_DISTANCE=  250.05

!
! SEGMENT 109
!
SEGMENT=  1483 1969 4453 4647
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00073 -0.00120
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00101 0.22429 -0.97452
SEGMENT_ORGX=  1233.24
SEGMENT_ORGY=  4550.45
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 110
!
SEGMENT=  1977 2463 4453 4647
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00052 -0.00035
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00022 0.22472 -0.97442
SEGMENT_ORGX=  1232.58
SEGMENT_ORGY=  4549.57
SEGMENT_DISTANCE=  250.06

!
! SEGMENT 111
!
SEGMENT=  1 487 4665 4859
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00029 0.00012
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00010 0.08035 -0.99677
SEGMENT_ORGX=  1231.02
SEGMENT_ORGY=  4762.17
SEGMENT_DISTANCE=  250.09

!
! SEGMENT 112
!
SEGMENT=  495 981 4665 4859
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00088 -0.00026
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00033 0.08081 -0.99673
SEGMENT_ORGX=  1230.03
SEGMENT_ORGY=  4761.35
SEGMENT_DISTANCE=  250.16

!
! SEGMENT 113
!
SEGMENT=  989 1475 4665 4859
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00118 -0.00039
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00029 0.07994 -0.99680
SEGMENT_ORGX=  1233.32
SEGMENT_ORGY=  4762.71
SEGMENT_DISTANCE=  250.02

!
! SEGMENT 114
!
SEGMENT=  1483 1969 4665 4859
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00068 0.00041
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00047 0.07984 -0.99681
SEGMENT_ORGX=  1232.24
SEGMENT_ORGY=  4762.45
SEGMENT_DISTANCE=  250.13

!
! SEGMENT 115
!
SEGMENT=  1977 2463 4665 4859
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00057 -0.00031
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00026 0.08104 -0.99671
SEGMENT_ORGX=  1232.41
SEGMENT_ORGY=  4761.16
SEGMENT_DISTANCE=  250.08

!
! SEGMENT 116
!
SEGMENT=  1 487 4877 5071
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00100 0.00013
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00007 -0.06472 -0.99790
SEGMENT_ORGX=  1232.66
SEGMENT_ORGY=  4973.59
SEGMENT_DISTANCE=  249.77

!
! SEGMENT 117
!
SEGMENT=  495 981 4877 5071
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00023 -0.00039
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00038 -0.06512 -0.99788
SEGMENT_ORGX=  1230.95
SEGMENT_ORGY=  4973.97
SEGMENT_DISTANCE=  250.03

!
! SEGMENT 118
!
SEGMENT=  989 1475 4877 5071
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00054 0.00036
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00033 -0.06446 -0.99792
SEGMENT_ORGX=  1232.19
SEGMENT_ORGY=  4973.06
SEGMENT_DISTANCE=  249.99

!
! SEGMENT 119
!
SEGMENT=  1483 1969 4877 5071
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 0.00061 -0.00046
DIRECTION_OF_SEGMENT_Y-AXIS=  -0.00041 -0.06543 -0.99786
SEGMENT_ORGX=  1230.82
SEGMENT_ORGY=  4974.83
SEGMENT_DISTANCE=  249.85

!
! SEGMENT 120
!
SEGMENT=  1977 2463 4877 5071
DIRECTION_OF_SEGMENT_X-AXIS=  1.00000 -0.00049 0.00007
DIRECTION_OF_SEGMENT_Y-AXIS=  0.00004 -0.06464 -0.99791
SEGMENT_ORGX=  1232.58
SEGMENT_ORGY=  4973.30
SEGMENT_DISTANCE=  249.98

"""
        )

except Exception:
    print("\n    ! Something went wrong when attempting to write XDS.INP . Aborting !\n")
    quit()

print("\n    !!!!!!!!!!!!!!!!!!!!!!!\n    !! Generated XDS.INP !! \n    !!!!!!!!!!!!!!!!!!!!!!! \n")

try:
    shutil.copy("XDS.INP", "XDS.INP.original")
except Exception:
    print("\n    ! Something went wrong when attempting to make a copy XDS.INP !\n")
    quit()
