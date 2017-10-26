# mitoseg
Mitochondria Segmentation Tool

Here are the command line options for the program:

    mitoseg -options <parameters> <filename pattern>

Options:

    -zrange <start slice #> <end slice #> Specify z-range. This options specifies z-range by two slice #. (the numbers that are contained in the filename of dataset slice images.

    -psize <pixel size>            Specify pixel size as nm/px.

    -roi <top> <left> <width> <height> Specify region of interest. Top, left pixel coordinates and width,height of bounding rectangle that you are interested in segmenting. If you define a roi, the program will work much faster.

    -src <source directory>            Specify source directory. Source directory is the path name of the dataset image files.

    -dst <destination directory>        Specify destination directory. This directory is used to store output files. The program creates many files regarding intermediate steps. So, I recommend to create a new "work" folder for segmentation purposes.

    -phase <phase #>            Apply phase <phase #> only. (must be 1, 2 or 3). The method consists of 3 phases. If you do not specify this option, all phases will be executed in order. The "-valid" and "-thick" options below affect 2nd and 3rd phases respectively. Hence, if you want to re-segment the same portion of the dataset by only changing these parameters, you may want to restart a particular phase manually in order to make it faster by avoiding redundant computation.

    -valid <validity threshold>        Specify validity threshold between 0-1 (default: 0.75). This is the validity threshold parameter explained in our last paper.

    -thick <z-thickness>            Specify snake z-thickness between 5-500 (default: 20), set 'full' to use the whole z-range specified in "-zrange" parameter. Snake thickness affects the segmentation accuracy. Thicker snakes provide continuous smooth segmentation along the z-range, but it increases false negatives. Default number is 20 and it generally produces better results w.r.t. our experiments.

    -cores <# of cpu cores>            Set # of cpu cores to utilize (default: 1). If your computer has multiple cpus or multiple cores, you can deploy some or all of them to work in collaboration in order to speed-up the process.

The options -zrange, -psize and <filename_pattern> are mandatory input.

Example:

    mitoseg -zrange 30 100 -psize 2.0 dataset_slice%04d.tif
    
        Process files dataset_slice0030.tif ... dataset_slice0100.tif assuming that pixel size is 2.0nm
        
    mitoseg -zrange 40 120 -psize 1.1 mito%d.bmp
    
        Process files mito40.bmp ... mito120.bmp assuming that pixel size is 1.1nm

Filename pattern: This pattern contains file name of slice images and slice number tag. "%d" indicates the position of the slice number tag in the filename. For example, "mito%d.bmp" generates file names such as mito0.bmp, mito1.bmp, mito12.bmp ... mito123.bmp and so on.
"mito%03d.bmp" generates file names by adding zeros to slice # to make it 3 digits. 000, 001, 002... and so on.

Many image types are supported (bmp, jpg, png, tiff...)

The program produces intermediate data and image files, final segmentation image files, Imod model file (.mod) and 3D mesh file (.ply) which can be visualized with an appropriate software such as Imod and MeshLab (free).
