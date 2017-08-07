#include "comm.h"
#include <iostream>
#include "fg_atoms.h"
#include "mapping.h"
#include "matrix_C.h"
#include <cstring>
#include <fstream>
#include "pointers.h"
#include "assert.h"


Fg_atoms::Fg_atoms(Mapping *map) : Pointers(map)
{
    input = new FrameSource();
}

void Fg_atoms::init(int argc, char **argv)
{
    fgtrj.open("FG_TRJ.lmpstrj", std::ofstream::out);
    fgtrj << "FG_TRJ : Coloring by C[5][i]" << std::endl;

    nframes = atoi(argv[2]);
    nSkipFrames = 0;
    if (argc > 3)
        nSkipFrames = atoi(argv[5]);

    printf("starting from Frame %d\n", nSkipFrames);

    assert(argv[1] != NULL);

    strcpy(input->trajName, argv[1]);

    input->trajectory.open(input->trajName, std::ifstream::in);
    if (input->trajectory.fail())
    {
        printf("Fail to load lammps trajectory file %s !\n", input->trajName);
        exit(-1);
    }

    readLammpsHeader();

    if (fg_num == 0)
    {
        printf("Error happened when parsing the trajectory, exit...\n");
        exit(-1);
    }

    r = new double* [fg_num];
    for (int i = 0; i < fg_num; i++)
    {
        r[i] = new double[3];
        r[i][0] = 0.0;
        r[i][1] = 0.0;
        r[i][2] = 0.0;
    }

    v = new double [3 * fg_num];
    memset(v, 0, sizeof(double) * 3 * fg_num);

    readInitialFrame();

    for (int i = 0; i < nSkipFrames; i++)
    {
        readNextFrame();
        printf("reading frame %d\r", i + 1);
    }
}

void Fg_atoms::readInitialFrame()
{
    input->elements = new std::string[input->header_size];
    readLammpsBody();
}

void Fg_atoms::readNextFrame()
{
    readLammpsHeader();
    readLammpsBody();
}

void Fg_atoms::cleanup()
{
    fgtrj.close();
    for (int i = 0; i < fg_num; i++) delete[] r[i];
    delete[] r;
    delete[] v;
}

void Fg_atoms::output()
{
    double **C = matrix_C->C;

    fgtrj << "ITEM: TIMESTEP" << std::endl;
    fgtrj << currentStep << std::endl;
    fgtrj << "ITEM: NUMBER OF ATOMS" << std::endl;
    fgtrj << fg_num << std::endl;
    fgtrj << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
    for (int i = 0; i < 3; i++) fgtrj << "0 " << L << std::endl;
    fgtrj << "ITEM: ATOMS id type m x y z vx vy vz" << std::endl;

    for (int i = 0; i < fg_num; i++)
    {
        fgtrj << i + 1 << ' ' << 1 << ' ' << 1 << ' ' << r[i][0] << ' ' << r[i][1] << ' ' << r[i][2] << ' ' << C[4][i] << ' ' << 0 << ' ' << 0 << std::endl;
    }
}

void Fg_atoms::readLammpsHeader()
{
    double low = 0.0;
    double high = 0.0;
    std::string line;
    int flag = 1;

    while (flag == 1)
    {
        // read next line in the header
        std::getline(input->trajectory, line);

        if (line.compare(0, 5, "ITEM:") == 0)
        {
            if (line.compare(6, 15, "NUMBER OF ATOMS") == 0)
            {
                input->trajectory >> fg_num;
            }
            else if (line.compare(6, 10, "BOX BOUNDS") == 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    input->trajectory >> low >> high;
                }
                L = high - low;
            }
            else if (line.compare(6, 8, "TIMESTEP") == 0)
            {
                input->trajectory >> currentStep;
            }
            else if (line.compare(6, 5, "ATOMS") == 0)
            {
                flag = 0;
                size_t prev = 11;
                size_t next = 0;
                int set_x = 0;
                int set_v = 0;
                int set_type = 0;
                input->header_size = 0;

                //Check if the input trajectory is valid or not
                while ((next = line.find_first_of(" ", prev)) != std::string::npos)
                {
                    if ((next - prev) == 0)
                    {
                        prev++;
                        continue;
                    }
                    else if (line.compare(prev, 1, "x") == 0)
                    {
                        input->x_pos = input->header_size;
                        set_x = 1;
                    }
                    else if (line.compare(prev, 2, "vx") == 0)
                    {
                        input->v_pos = input->header_size;
                        set_v = 1;
                    }
                    else if (line.compare(prev, 4, "type") == 0)
                    {
                        input->type_pos = input->header_size;
                        set_type = 1;
                    }
                    input->header_size++;
                    prev = next;
                }

                if (prev < line.size())
                {
                    if (line.compare(prev, 1, "x") == 0)
                    {
                        input->x_pos = input->header_size;
                        set_x = 1;
                    }
                    else if (line.compare(prev, 2, "vx") == 0)
                    {
                        input->v_pos = input->header_size;
                        set_v = 1;
                    }
                    else if (line.compare(prev, 4, "type") == 0)
                    {
                        input->type_pos = input->header_size;
                        set_type = 1;
                    }
                    input->header_size++;
                }

                if ((set_x == 0) || (set_v == 0))
                {
                    printf("Warning: Either position of speed was not set\n");
                    exit(-1);
                }
            }
            else
            {
                printf("Unrecognized line inheader %s", line.c_str());
            }
        }
    }
    return;
}

void Fg_atoms::readLammpsBody()
{
    int j = 0;
    std::string line;

    for (int i = 0; i < fg_num; i++) // totally fg_num lines in the file
    {
        std::getline(input->trajectory, line, '\n');
        if ((j = stringSplit(line, " \t", input->elements)) != input->header_size)
        {
            printf("Warning: Number of fields detedted in frame body");
            printf("(%d) does not agree with number expected from frame header (%d)!\n", j, input->header_size);
        }

        for (j = 0; j < 3; j++)
        {
            r[i][j] = atof(input->elements[j + input->x_pos].c_str());
        }

        for (j = 0; j < 3; j++)
        {
            v[i + j * fg_num] = atof(input->elements[j + input->v_pos].c_str());
        }

    }
}

int Fg_atoms::stringSplit(std::string source, const char *const delimiter, std::string *results)
{
    int count = 0;
    size_t prev = 0;
    size_t next = 0;

    while ((next = source.find_first_of(delimiter, prev)) != std::string::npos)
    {
        if (next - prev != 0) //not the end
        {
            results[count] = source.substr(prev, next - prev);
            count++;
        }
        prev = next + 1;
    }

    if (prev < source.size())  //something left, the last elements may end up without space or tab
    {
        results[count] = source.substr(prev);
        count++;
    }
    return count;
}

void Fg_atoms::finishReading()
{
    input->trajectory.close();
    delete[] input->elements;
}

