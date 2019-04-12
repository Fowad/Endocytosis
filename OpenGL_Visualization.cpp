#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <assert.h>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Shader.h"

using namespace std;

const int n_filaments = 144, n_central = 36; // BEWARE OF THE LINE N_MAX
//const int n_gel = 25; // n_gel = 2*n_filaments + 1; number of filament bars plus number of gaps between them
const int m = 12, m_p = 6;

double f_stall = 7.3; //pN
double f_max_each_puller = ((n_filaments - n_central)*f_stall)/n_central;

const unsigned long long int Ntstep = 7*pow(10,9); // number of time steps.
const unsigned long long int Ntstep_to_store_data = pow(10,6);

double opengl_window_scaling_metric = 50; // meaning 1 unit in the opengl window equals that number in nanometer. And the point (0,0) is at the center!

int main( )
{
    // Window dimensions
    const GLuint WIDTH = 600, HEIGHT = 600;
    
    // Init GLFW
    glfwInit();
    
    // Set all the required options for GLFW
    glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 4 );
    glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 1 );
    glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
    glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE );
    glfwWindowHint( GLFW_RESIZABLE, GL_FALSE );
    
    // Create a GLFWwindow object that we can use for GLFW's functions
    GLFWwindow *window = glfwCreateWindow( WIDTH, HEIGHT, "Result", nullptr, nullptr );
    
    int screenWidth, screenHeight;
    glfwGetFramebufferSize( window, &screenWidth, &screenHeight );
    
    if ( nullptr == window )
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate( );
        return EXIT_FAILURE;
    }
    
    glfwMakeContextCurrent( window );
    
    // Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    glewExperimental = GL_TRUE;
    
    // Initialize GLEW to setup the OpenGL Function pointers
    if ( GLEW_OK != glewInit( ) )
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }
    
    glViewport( 0, 0, screenWidth, screenHeight ); // Define the viewport dimensions
    
    Shader ourShader( "res/shaders/core.vs", "res/shaders/core.frag" ); // Build and compile our shader program

    GLdouble Force_array [m][m][Ntstep/Ntstep_to_store_data];
    //GLdouble Force_array_final [m][m];
    
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            ifstream Force_data ("res/Force_vs_time_" + to_string(row) + "_" + to_string(col) + ".txt");

            if (Force_data.is_open()) {for (int i=0; i<Ntstep/Ntstep_to_store_data; i++){Force_data >> Force_array [row][col][i];}}}}
    

//        for (int row = 0; row < m; row++)
//        {
//            for (int col = 0; col < m; col++)
//            {
//                ifstream Force_dist ("res/Third_Column.txt");
//
//                if (Force_dist.is_open()) {Force_dist >> Force_array_final [row][col];}}}
    
    
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            for (int i=0; i<Ntstep/Ntstep_to_store_data; i++){

            if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) ){

                Force_array [row][col][i] = -1*Force_array [row][col][i]/f_max_each_puller;

            }

            else {

                Force_array [row][col][i] = Force_array [row][col][i]/f_stall;
            }

            //cout << Force_array_final [row][col] << endl;
                
        }}}

    
    GLuint VBO_for_a_point[1], VAO_for_a_point[1]; // the purpose of this point is to make a delay in drawing objects in the output window!
    glGenVertexArrays (1, VAO_for_a_point);
    glGenBuffers(1, VBO_for_a_point);
    
    GLuint VBOs_for_the_filaments [n_filaments];
    GLuint VAOs_for_the_filaments [n_filaments];
    GLuint EBOs_for_the_filaments [n_filaments];
    
    glGenVertexArrays(n_filaments, VAOs_for_the_filaments);
    glGenBuffers(n_filaments, VBOs_for_the_filaments);
    glGenBuffers(n_filaments, EBOs_for_the_filaments);

    const GLfloat n_max = 12; // maximum number of filaments that I would possibly draw in this opengl window.
    GLfloat delta_prime = 8.5/opengl_window_scaling_metric;
    GLfloat R = (2 - n_max*delta_prime)/(n_max + 1);
    GLfloat L = (m*delta_prime + (m-1)*R);
    GLfloat Sla2_Radius = 0.5*(n_central*delta_prime + (n_central+1)*R);
    
    /////////////////////////////////////////////
    //============== Game Loop ================//
    /////////////////////////////////////////////
    
    while ( !glfwWindowShouldClose( window ) )
    {
        for (int i=0; i<(Ntstep/Ntstep_to_store_data); i+=100){

        glfwPollEvents( ); // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
            
        // Render
        // Clear the korbuffer
        glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
        glClear( GL_COLOR_BUFFER_BIT );

//      Usually when you have multiple objects you want to draw,
//      you first generate/configure all the VAOs (and thus the required VBO and attribute pointers)
//      and store those for later use. The moment we want to draw one of our objects,
//      we take the corresponding VAO, bind it, then draw the object and unbind the VAO again.
            
        GLfloat filaments_vertices [m][m][24];
        GLuint filaments_indices [m][m][6];
            
        for (int row = 0; row < m; row++){for (int col = 0; col < m; col++){

            if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) ){

                filaments_vertices [row][col][0] = -0.5*L + (float(col)*(delta_prime+R)) + delta_prime; // Top right x
                filaments_vertices [row][col][1] = -0.5*L + (float(row)*(delta_prime+R)) + delta_prime; // Top right y
                filaments_vertices [row][col][2] = 0.0f;
                filaments_vertices [row][col][3] = float(Force_array[row][col][i]);
                filaments_vertices [row][col][4] = 0.0f;
                filaments_vertices [row][col][5] = 0.0f;
                filaments_vertices [row][col][6] = -0.5*L + (float(col)*(delta_prime+R)) + delta_prime; // Bottom right x
                filaments_vertices [row][col][7] = -0.5*L + (float(row)*(delta_prime+R)); // Bottom right y
                filaments_vertices [row][col][8] = 0.0f;
                filaments_vertices [row][col][9] = float(Force_array[row][col][i]);
                filaments_vertices [row][col][10] = 0.0f;
                filaments_vertices [row][col][11] = 0.0f;
                filaments_vertices [row][col][12] = -0.5*L + (float(col)*(delta_prime+R)); // Bottom Left x
                filaments_vertices [row][col][13] = -0.5*L + (float(row)*(delta_prime+R)); // Bottom left y
                filaments_vertices [row][col][14] = 0.0f;
                filaments_vertices [row][col][15] = float(Force_array[row][col][i]);
                filaments_vertices [row][col][16] = 0.0f;
                filaments_vertices [row][col][17] = 0.0f;
                filaments_vertices [row][col][18] = -0.5*L + (float(col)*(delta_prime+R)); // Top left x
                filaments_vertices [row][col][19] = -0.5*L + (float(row)*(delta_prime+R)) + delta_prime; // Top left y
                filaments_vertices [row][col][20] = 0.0f;
                filaments_vertices [row][col][21] = float(Force_array[row][col][i]);
                filaments_vertices [row][col][22] = 0.0f;
                filaments_vertices [row][col][23] = 0.0f;

                }

                else {
                    filaments_vertices [row][col][0] = -0.5*L + (float(col)*(delta_prime+R)) + delta_prime; // Top right x
                    filaments_vertices [row][col][1] = -0.5*L + (float(row)*(delta_prime+R)) + delta_prime; // Top right y
                    filaments_vertices [row][col][2] = 0.0f;
                    filaments_vertices [row][col][3] = 0.0f;
                    filaments_vertices [row][col][4] = 0.0f;
                    filaments_vertices [row][col][5] = float(Force_array[row][col][i]);
                    filaments_vertices [row][col][6] = -0.5*L + (float(col)*(delta_prime+R)) + delta_prime; // Bottom right x
                    filaments_vertices [row][col][7] = -0.5*L + (float(row)*(delta_prime+R)); // Bottom right y
                    filaments_vertices [row][col][8] = 0.0f;
                    filaments_vertices [row][col][9] = 0.0f;
                    filaments_vertices [row][col][10] = 0.0f;
                    filaments_vertices [row][col][11] = float(Force_array[row][col][i]);
                    filaments_vertices [row][col][12] = -0.5*L + (float(col)*(delta_prime+R)); // Bottom Left x
                    filaments_vertices [row][col][13] = -0.5*L + (float(row)*(delta_prime+R)); // Bottom left y
                    filaments_vertices [row][col][14] = 0.0f;
                    filaments_vertices [row][col][15] = 0.0f;
                    filaments_vertices [row][col][16] = 0.0f;
                    filaments_vertices [row][col][17] = float(Force_array[row][col][i]);
                    filaments_vertices [row][col][18] = -0.5*L + (float(col)*(delta_prime+R)); // Top left x
                    filaments_vertices [row][col][19] = -0.5*L + (float(row)*(delta_prime+R)) + delta_prime; // Top left y
                    filaments_vertices [row][col][20] = 0.0f;
                    filaments_vertices [row][col][21] = 0.0f;
                    filaments_vertices [row][col][22] = 0.0f;
                    filaments_vertices [row][col][23] = float(Force_array[row][col][i]);
                    }
        }}


            for (int row = 0; row < m; row++){for (int col = 0; col < m; col++){
                filaments_indices [row][col][0] = 0;
                filaments_indices [row][col][1] = 1;
                filaments_indices [row][col][2] = 2;
                filaments_indices [row][col][3] = 2;
                filaments_indices [row][col][4] = 3;
                filaments_indices [row][col][5] = 0;}}   
    
        // The fourth parameter specifies how we want the graphics card to manage the given data. This can take 3 forms:
        //GL_STATIC_DRAW: the data will most likely not change at all or very rarely.
        //GL_DYNAMIC_DRAW: the data is likely to change a lot.
        //GL_STREAM_DRAW: the data will change every time it is drawn.

            
        GLfloat test_filaments [24];
        GLuint test_i_filaments [6];
            
        for (int row = 0; row < m; row++){for (int col = 0; col < m; col++){
            
            for (int j=0; j<24; j++){test_filaments [j] = filaments_vertices[row][col][j];}
            for (int j=0; j<6; j++){test_i_filaments [j] = filaments_indices[row][col][j];}
        
            glBindVertexArray(VAOs_for_the_filaments[col+m*row]);
            glBindBuffer(GL_ARRAY_BUFFER, VBOs_for_the_filaments[col+m*row]);
            glBufferData(GL_ARRAY_BUFFER, sizeof(test_filaments), test_filaments, GL_STATIC_DRAW);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOs_for_the_filaments[col+m*row]);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(test_i_filaments), test_i_filaments, GL_STATIC_DRAW);
            glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * ) 0 );
            glEnableVertexAttribArray( 0 );
            glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
            glEnableVertexAttribArray( 1 );
            glBindVertexArray( 0 );}}

            ourShader.Use();

        for (int row = 0; row < m; row++){for (int col = 0; col < m; col++){
            glBindVertexArray(VAOs_for_the_filaments[col+m*row]);
            glDrawElements( GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0 );
            glBindVertexArray(0);}}

            glfwSwapBuffers( window ); // Swap the screen buffers

            for (int j=0; j<=500000; j++){

                if ( j < 500000 )
                {
                GLfloat point[] = { 0.0f, 0.0f, 0.0f,  0.5f, 0.5f, 0.0f };

                glBindVertexArray(VAO_for_a_point[0]);
                glBindBuffer(GL_ARRAY_BUFFER, VBO_for_a_point[0]);
                glBufferData(GL_ARRAY_BUFFER, sizeof(point), point, GL_STATIC_DRAW);
                glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * ) 0 );
                glEnableVertexAttribArray( 0 );
                glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
                glEnableVertexAttribArray( 1 );
                glBindVertexArray( 0 );

                glBindVertexArray(VAO_for_a_point[0]);
                glDrawArrays(GL_POINT, 0, 3);
                glBindVertexArray(0);
                }}

            cout << "t = " << i*(0.001/2) << " sec" << endl;
            
        }
    }
    
    // Properly de-allocate all resources once they've outlived their purpose
    glDeleteBuffers(1, VBO_for_a_point);
    glDeleteBuffers(n_filaments, VBOs_for_the_filaments);

    glDeleteVertexArrays(1, VAO_for_a_point);
    glDeleteVertexArrays(n_filaments, VAOs_for_the_filaments);

    glDeleteBuffers(n_filaments, EBOs_for_the_filaments);

    // Terminate GLFW, clearing any resources allocated by GLFW.
    glfwTerminate( );
    
    return EXIT_SUCCESS;
}
