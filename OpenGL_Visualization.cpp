#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <assert.h>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Camera.h"
#include "Shader.h"

using namespace std;

const int n_filaments = 64, n_central = 16; // BEWARE OF THE LINE N_MAX

const int m = 8, m_p = 4; // m = sqrt(n_filaments); m_p = sqrt(n_central)

const unsigned long long int Ntstep = 2*pow(10,9);
const unsigned long long int Ntstep_to_store_data = pow(10,6);

GLfloat opengl_window_scaling_metric = 100; // meaning 1 unit in the opengl window equals that number in nanometer. And the point (0,0) is at the center!

const GLuint WIDTH = 800, HEIGHT = 700;

int SCREEN_WIDTH, SCREEN_HEIGHT;

void KeyCallback ( GLFWwindow *window, int key, int scancode, int action, int mode );
void ScrollCallback ( GLFWwindow *window, double xOffset, double yOffset );
void MouseCallback ( GLFWwindow *window, double xPos, double yPos );
void DoMovement();

Camera camera( glm::vec3 (0.0f, 0.0f, 3.0f) );
GLfloat lastX = WIDTH / 2.0f;
GLfloat lastY = WIDTH / 2.0f;
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

bool keys[1024];
bool firstMouse = true;

int main( )
{
    
    glfwInit();
    glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 4 );
    glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 1 );
    glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
    glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE );
    glfwWindowHint( GLFW_RESIZABLE, GL_FALSE );
    
    GLFWwindow *window = glfwCreateWindow( WIDTH, HEIGHT, "Result", nullptr, nullptr );
    
    glfwGetFramebufferSize( window, &SCREEN_WIDTH, &SCREEN_HEIGHT );
    
    if ( nullptr == window )
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate( );
        return EXIT_FAILURE;
    }
    
    glfwMakeContextCurrent( window );
    glfwSetKeyCallback( window, KeyCallback );
    glfwSetCursorPosCallback( window, MouseCallback );
    glfwSetScrollCallback( window, ScrollCallback );
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glewExperimental = GL_TRUE;
    
    if ( GLEW_OK != glewInit( ) )
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }
    
    glViewport( 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT ); // Define the viewport dimensions
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    
    
    Shader ourShader( "res/shaders/core.vs", "res/shaders/core.frag" ); // Build and compile our shader program

    /////////////////////////////////////////////
    //============== Arrays ================////
    /////////////////////////////////////////////
    
    const GLfloat n_max = 10; // maximum number of filaments that I would possibly draw in this opengl window.
    GLfloat delta_prime = 2.7/opengl_window_scaling_metric;
    GLfloat d = delta_prime/2.0f;
    GLfloat R = (2 - n_max*delta_prime)/(n_max + 1);
    GLfloat L = (m*delta_prime + (m-1)*R);
//    GLfloat Sla2_Radius = 0.5*(n_central*delta_prime + (n_central+1)*R);
    
    GLdouble z_sc [Ntstep/Ntstep_to_store_data];

    GLdouble Num_subunits_array [m][m][Ntstep/Ntstep_to_store_data], dz_bend_array [m][m][Ntstep/Ntstep_to_store_data], dz_elas_array [m][m][Ntstep/Ntstep_to_store_data];
    
    ifstream Sc_zdata_prime ("res/Scaled_zdata.txt");
    if (Sc_zdata_prime.is_open()){
    for (int i=0;i<Ntstep/Ntstep_to_store_data;i++){Sc_zdata_prime >> z_sc[i];}}

    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){ifstream Num_subunits_prime ("res/Num_subunits_" + to_string(r) + "_" + to_string(c) + ".txt"), dz_bend_prime ("res/dz_bend_" + to_string(r) + "_" + to_string(c) + ".txt"), dz_elas_prime ("res/dz_elas_" + to_string(r) + "_" + to_string(c) + ".txt"); if (Num_subunits_prime.is_open() && dz_bend_prime.is_open() && dz_elas_prime.is_open()) {for (int i=0; i<Ntstep/Ntstep_to_store_data; i++){Num_subunits_prime >> Num_subunits_array[r][c][i]; dz_bend_prime >> dz_bend_array[r][c][i]; dz_elas_prime >> dz_elas_array[r][c][i];}}}}

    
    GLuint VBO_for_the_mem, VAO_for_the_mem;
    glGenVertexArrays (1, &VAO_for_the_mem);
    glGenBuffers(1, &VBO_for_the_mem);
    
    GLuint VBO_for_a_point, VAO_for_a_point; // the purpose of this point is to make a delay in drawing objects in the final window
    glGenVertexArrays (1, &VAO_for_a_point);
    glGenBuffers(1, &VBO_for_a_point);
    
    GLuint VBOs_for_the_filaments [n_filaments];
    GLuint VAOs_for_the_filaments [n_filaments];
    
    glGenVertexArrays(n_filaments, VAOs_for_the_filaments);
    glGenBuffers(n_filaments, VBOs_for_the_filaments);
    
    // Cube examples
    GLuint VBO, VAO;
    glGenVertexArrays( 1, &VAO );
    glGenBuffers( 1, &VBO );
    
    
//
//    /////////////////////////////////////////////
//    //============== Game Loop ================//
//    /////////////////////////////////////////////
    
    
    
    while ( !glfwWindowShouldClose( window ) ){
        
        for (int i=0; i<(Ntstep/Ntstep_to_store_data); i++){
        
        // Set frame time
        GLfloat currentFrame = glfwGetTime( );
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        glfwPollEvents( );
        DoMovement();
            
        glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        
            
//      Usually when you have multiple objects you want to draw,
//      you first generate/configure all the VAOs (and thus the required VBO and attribute pointers)
//      and store those for later use. The moment we want to draw one of our objects,
//      we take the corresponding VAO, bind it, then draw the object and unbind the VAO again.
            
            
            
//        GLfloat membrane_vertices [] =
//        {
//            // Positions          // Colors
//            0.95f, -0.75f + (3.0f/5.0f)*delta_prime + float(z_sc[i]), 0.0f,    0.0f, 0.0f, 1.0f, // Top Right
//            0.95f, -0.75f + float(z_sc[i]), 0.0f,   0.0f, 0.0f, 1.0f, // Bottom Right
//            -0.95f, -0.75f + float(z_sc[i]), 0.0f,   0.0f, 0.0f, 1.0f, // Bottom Left
//            -0.95f, -0.75f + (3.0f/5.0f)*delta_prime + float(z_sc[i]), 0.0f,   0.0f, 0.0f, 1.0f // Top Left
//        };
//        GLuint membrane_indices [] =
//        {
//                0, 1, 2, // First Triangle
//                2, 3, 0  // Second Triangle
//        };
            
            
//            glBindVertexArray(VAO_for_the_mem[0]);
//            glBindBuffer(GL_ARRAY_BUFFER, VBO_for_the_mem[0]);
//            glBufferData(GL_ARRAY_BUFFER, sizeof(membrane_vertices), membrane_vertices, GL_STREAM_DRAW);
//            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_for_the_mem[0]);
//            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(membrane_indices), membrane_indices, GL_STREAM_DRAW);
//            // Position attribute
//            glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * ) 0 );
//            glEnableVertexAttribArray( 0 );
//            // Color attribute
//            glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
//            glEnableVertexAttribArray( 1 );
//            glBindVertexArray( 0 ); // Unbind VAO
            
            
            
            
            GLfloat membrane_vertices [] =
            {
                
                -1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f, // front
                1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f  + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                
                -1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f, // back
                1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f, // top
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                
                -1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,// bottom
                1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                
                1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f, // right
                1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,
                
                -1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f,// left
                -1.0f, -1.0f + float(z_sc[i]), -1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, -1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]) + (3.0f/5.0f)*delta_prime, 1.0f,  0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f + float(z_sc[i]), 1.0f,  0.0f, 0.0f, 1.0f
                
                
            };

            
            // Bind our Vertex Array Object first, then bind and set our buffers and pointers.
            glBindVertexArray( VAO_for_the_mem );
            
            glBindBuffer( GL_ARRAY_BUFFER, VBO_for_the_mem );
            glBufferData(GL_ARRAY_BUFFER, sizeof(membrane_vertices), membrane_vertices, GL_STREAM_DRAW);
            
            // Position attribute
            glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )0 );
            glEnableVertexAttribArray( 0 );
            // Color attribute
            glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
            glEnableVertexAttribArray( 1 );
            
            glBindVertexArray( 0 ); // Unbind VAO

            
            
          
            GLfloat filaments_vertices [m][m][216]; // Pay attention that for the triangles in the top and bottom parts, the order of 1 2 3 is from top left to top right to bottom right. but for the other triangles the order is bottom left to bottom right to top right.

            
        for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
            
            if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) )
            
            {
            
                filaments_vertices [r][c][0] = -0.5*L + (float(c)*(delta_prime+R)); // front
                filaments_vertices [r][c][1] = -1.0f;
                filaments_vertices [r][c][2] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][3] = 0.0f;
                filaments_vertices [r][c][4] = 1.0f;
                filaments_vertices [r][c][5] = 0.0f;
                filaments_vertices [r][c][6] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][7] = -1.0f;
                filaments_vertices [r][c][8] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][9] = 0.0f;
                filaments_vertices [r][c][10] = 1.0f;
                filaments_vertices [r][c][11] = 0.0f;
                filaments_vertices [r][c][12] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][13] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][14] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][15] = 0.0f;
                filaments_vertices [r][c][16] = 1.0f;
                filaments_vertices [r][c][17] = 0.0f;
                filaments_vertices [r][c][18] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][19] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][20] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][21] = 0.0f;
                filaments_vertices [r][c][22] = 1.0f;
                filaments_vertices [r][c][23] = 0.0f;
                filaments_vertices [r][c][24] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][25] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][26] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][27] = 0.0f;
                filaments_vertices [r][c][28] = 1.0f;
                filaments_vertices [r][c][29] = 0.0f;
                filaments_vertices [r][c][30] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][31] = -1.0f;
                filaments_vertices [r][c][32] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][33] = 0.0f;
                filaments_vertices [r][c][34] = 1.0f;
                filaments_vertices [r][c][35] = 0.0f;
                
                filaments_vertices [r][c][36] = -0.5*L + (float(c)*(delta_prime+R)); // back
                filaments_vertices [r][c][37] = -1.0f;
                filaments_vertices [r][c][38] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][39] = 0.0f;
                filaments_vertices [r][c][40] = 1.0f;
                filaments_vertices [r][c][41] = 0.0f;
                filaments_vertices [r][c][42] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][43] = -1.0f;
                filaments_vertices [r][c][44] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][45] = 0.0f;
                filaments_vertices [r][c][46] = 1.0f;
                filaments_vertices [r][c][47] = 0.0f;
                filaments_vertices [r][c][48] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][49] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][50] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][51] = 0.0f;
                filaments_vertices [r][c][52] = 1.0f;
                filaments_vertices [r][c][53] = 0.0f;
                filaments_vertices [r][c][54] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][55] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][56] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][57] = 0.0f;
                filaments_vertices [r][c][58] = 1.0f;
                filaments_vertices [r][c][59] = 0.0f;
                filaments_vertices [r][c][60] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][61] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][62] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][63] = 0.0f;
                filaments_vertices [r][c][64] = 1.0f;
                filaments_vertices [r][c][65] = 0.0f;
                filaments_vertices [r][c][66] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][67] = -1.0f;
                filaments_vertices [r][c][68] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][69] = 0.0f;
                filaments_vertices [r][c][70] = 1.0f;
                filaments_vertices [r][c][71] = 0.0f;
                
                filaments_vertices [r][c][72] = -0.5*L + (float(c)*(delta_prime+R)); // top
                filaments_vertices [r][c][73] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][74] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][75] = 0.0f;
                filaments_vertices [r][c][76] = 1.0f;
                filaments_vertices [r][c][77] = 0.0f;
                filaments_vertices [r][c][78] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][79] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][80] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][81] = 0.0f;
                filaments_vertices [r][c][82] = 1.0f;
                filaments_vertices [r][c][83] = 0.0f;
                filaments_vertices [r][c][84] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][85] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][86] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][87] = 0.0f;
                filaments_vertices [r][c][88] = 1.0f;
                filaments_vertices [r][c][89] = 0.0f;
                filaments_vertices [r][c][90] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][91] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][92] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][93] = 0.0f;
                filaments_vertices [r][c][94] = 1.0f;
                filaments_vertices [r][c][95] = 0.0f;
                filaments_vertices [r][c][96] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][97] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][98] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][99] = 0.0f;
                filaments_vertices [r][c][100] = 1.0f;
                filaments_vertices [r][c][101] = 0.0f;
                filaments_vertices [r][c][102] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][103] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][104] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][105] = 0.0f;
                filaments_vertices [r][c][106] = 1.0f;
                filaments_vertices [r][c][107] = 0.0f;
                
                filaments_vertices [r][c][108] = -0.5*L + (float(c)*(delta_prime+R)); // bottom
                filaments_vertices [r][c][109] = -1.0f;
                filaments_vertices [r][c][110] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][111] = 0.0f;
                filaments_vertices [r][c][112] = 1.0f;
                filaments_vertices [r][c][113] = 0.0f;
                filaments_vertices [r][c][114] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][115] = -1.0f;
                filaments_vertices [r][c][116] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][117] = 0.0f;
                filaments_vertices [r][c][118] = 1.0f;
                filaments_vertices [r][c][119] = 0.0f;
                filaments_vertices [r][c][120] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][121] = -1.0f;
                filaments_vertices [r][c][122] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][123] = 0.0f;
                filaments_vertices [r][c][124] = 1.0f;
                filaments_vertices [r][c][125] = 0.0f;
                filaments_vertices [r][c][126] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][127] = -1.0f;
                filaments_vertices [r][c][128] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][129] = 0.0f;
                filaments_vertices [r][c][130] = 1.0f;
                filaments_vertices [r][c][131] = 0.0f;
                filaments_vertices [r][c][132] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][133] = -1.0f;
                filaments_vertices [r][c][134] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][135] = 0.0f;
                filaments_vertices [r][c][136] = 1.0f;
                filaments_vertices [r][c][137] = 0.0f;
                filaments_vertices [r][c][138] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][139] = -1.0f;
                filaments_vertices [r][c][140] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][141] = 0.0f;
                filaments_vertices [r][c][142] = 1.0f;
                filaments_vertices [r][c][143] = 0.0f;
                
                filaments_vertices [r][c][144] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime; // right
                filaments_vertices [r][c][145] = -1.0f;
                filaments_vertices [r][c][146] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][147] = 0.0f;
                filaments_vertices [r][c][148] = 1.0f;
                filaments_vertices [r][c][149] = 0.0f;
                filaments_vertices [r][c][150] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][151] = -1.0f;
                filaments_vertices [r][c][152] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][153] = 0.0f;
                filaments_vertices [r][c][154] = 1.0f;
                filaments_vertices [r][c][155] = 0.0f;
                filaments_vertices [r][c][156] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][157] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][158] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][159] = 0.0f;
                filaments_vertices [r][c][160] = 1.0f;
                filaments_vertices [r][c][161] = 0.0f;
                filaments_vertices [r][c][162] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][163] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][164] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][165] = 0.0f;
                filaments_vertices [r][c][166] = 1.0f;
                filaments_vertices [r][c][167] = 0.0f;
                filaments_vertices [r][c][168] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][169] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][170] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][171] = 0.0f;
                filaments_vertices [r][c][172] = 1.0f;
                filaments_vertices [r][c][173] = 0.0f;
                filaments_vertices [r][c][174] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][175] = -1.0f;
                filaments_vertices [r][c][176] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][177] = 0.0f;
                filaments_vertices [r][c][178] = 1.0f;
                filaments_vertices [r][c][179] = 0.0f;
                
                filaments_vertices [r][c][180] = -0.5*L + (float(c)*(delta_prime+R)); // left
                filaments_vertices [r][c][181] = -1.0f;
                filaments_vertices [r][c][182] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][183] = 0.0f;
                filaments_vertices [r][c][184] = 1.0f;
                filaments_vertices [r][c][185] = 0.0f;
                filaments_vertices [r][c][186] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][187] = -1.0f;
                filaments_vertices [r][c][188] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][189] = 0.0f;
                filaments_vertices [r][c][190] = 1.0f;
                filaments_vertices [r][c][191] = 0.0f;
                filaments_vertices [r][c][192] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][193] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][194] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][195] = 0.0f;
                filaments_vertices [r][c][196] = 1.0f;
                filaments_vertices [r][c][197] = 0.0f;
                filaments_vertices [r][c][198] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][199] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][200] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][201] = 0.0f;
                filaments_vertices [r][c][202] = 1.0f;
                filaments_vertices [r][c][203] = 0.0f;
                filaments_vertices [r][c][204] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][205] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][206] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][207] = 0.0f;
                filaments_vertices [r][c][208] = 1.0f;
                filaments_vertices [r][c][209] = 0.0f;
                filaments_vertices [r][c][210] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][211] = -1.0f;
                filaments_vertices [r][c][212] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][213] = 0.0f;
                filaments_vertices [r][c][214] = 1.0f;
                filaments_vertices [r][c][215] = 0.0f;
                
                
                
               /* filaments_vertices [r][c][] =
                
                {
                    
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f, // front
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    
                    
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f, // back
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    
                    
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f, // top
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    
                    
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f, // bottom
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - delta_prime - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    
                    
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f, // right
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)) + delta_prime, -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f
                    
                    
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f, // left
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)) - delta_prime,  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f,
                    -0.5*L + (float(c)*(delta_prime+R)), -1.0f, 1.0f - (float(r)*(delta_prime+R)),  0.0f, 1.0f, 0.0f
                }; */
            
            }
            
            else {
            
                                
                filaments_vertices [r][c][0] = -0.5*L + (float(c)*(delta_prime+R)); // front
                filaments_vertices [r][c][1] = -1.0f;
                filaments_vertices [r][c][2] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][3] = 1.0f;
                filaments_vertices [r][c][4] = 0.0f;
                filaments_vertices [r][c][5] = 0.0f;
                filaments_vertices [r][c][6] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][7] = -1.0f;
                filaments_vertices [r][c][8] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][9] = 1.0f;
                filaments_vertices [r][c][10] = 0.0f;
                filaments_vertices [r][c][11] = 0.0f;
                filaments_vertices [r][c][12] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][13] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][14] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][15] = 1.0f;
                filaments_vertices [r][c][16] = 0.0f;
                filaments_vertices [r][c][17] = 0.0f;
                filaments_vertices [r][c][18] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][19] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][20] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][21] = 1.0f;
                filaments_vertices [r][c][22] = 0.0f;
                filaments_vertices [r][c][23] = 0.0f;
                filaments_vertices [r][c][24] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][25] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][26] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][27] = 1.0f;
                filaments_vertices [r][c][28] = 0.0f;
                filaments_vertices [r][c][29] = 0.0f;
                filaments_vertices [r][c][30] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][31] = -1.0f;
                filaments_vertices [r][c][32] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][33] = 1.0f;
                filaments_vertices [r][c][34] = 0.0f;
                filaments_vertices [r][c][35] = 0.0f;
                
                filaments_vertices [r][c][36] = -0.5*L + (float(c)*(delta_prime+R)); // back
                filaments_vertices [r][c][37] = -1.0f;
                filaments_vertices [r][c][38] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][39] = 1.0f;
                filaments_vertices [r][c][40] = 0.0f;
                filaments_vertices [r][c][41] = 0.0f;
                filaments_vertices [r][c][42] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][43] = -1.0f;
                filaments_vertices [r][c][44] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][45] = 1.0f;
                filaments_vertices [r][c][46] = 0.0f;
                filaments_vertices [r][c][47] = 0.0f;
                filaments_vertices [r][c][48] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][49] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][50] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][51] = 1.0f;
                filaments_vertices [r][c][52] = 0.0f;
                filaments_vertices [r][c][53] = 0.0f;
                filaments_vertices [r][c][54] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][55] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][56] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][57] = 1.0f;
                filaments_vertices [r][c][58] = 0.0f;
                filaments_vertices [r][c][59] = 0.0f;
                filaments_vertices [r][c][60] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][61] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][62] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][63] = 1.0f;
                filaments_vertices [r][c][64] = 0.0f;
                filaments_vertices [r][c][65] = 0.0f;
                filaments_vertices [r][c][66] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][67] = -1.0f;
                filaments_vertices [r][c][68] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][69] = 1.0f;
                filaments_vertices [r][c][70] = 0.0f;
                filaments_vertices [r][c][71] = 0.0f;
                
                filaments_vertices [r][c][72] = -0.5*L + (float(c)*(delta_prime+R)); // top
                filaments_vertices [r][c][73] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][74] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][75] = 1.0f;
                filaments_vertices [r][c][76] = 0.0f;
                filaments_vertices [r][c][77] = 0.0f;
                filaments_vertices [r][c][78] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][79] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][80] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][81] = 1.0f;
                filaments_vertices [r][c][82] = 0.0f;
                filaments_vertices [r][c][83] = 0.0f;
                filaments_vertices [r][c][84] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][85] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][86] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][87] = 1.0f;
                filaments_vertices [r][c][88] = 0.0f;
                filaments_vertices [r][c][89] = 0.0f;
                filaments_vertices [r][c][90] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][91] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][92] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][93] = 1.0f;
                filaments_vertices [r][c][94] = 0.0f;
                filaments_vertices [r][c][95] = 0.0f;
                filaments_vertices [r][c][96] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][97] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][98] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][99] = 1.0f;
                filaments_vertices [r][c][100] = 0.0f;
                filaments_vertices [r][c][101] = 0.0f;
                filaments_vertices [r][c][102] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][103] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][104] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][105] = 1.0f;
                filaments_vertices [r][c][106] = 0.0f;
                filaments_vertices [r][c][107] = 0.0f;
                
                filaments_vertices [r][c][108] = -0.5*L + (float(c)*(delta_prime+R)); // bottom
                filaments_vertices [r][c][109] = -1.0f;
                filaments_vertices [r][c][110] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][111] = 1.0f;
                filaments_vertices [r][c][112] = 0.0f;
                filaments_vertices [r][c][113] = 0.0f;
                filaments_vertices [r][c][114] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][115] = -1.0f;
                filaments_vertices [r][c][116] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][117] = 1.0f;
                filaments_vertices [r][c][118] = 0.0f;
                filaments_vertices [r][c][119] = 0.0f;
                filaments_vertices [r][c][120] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][121] = -1.0f;
                filaments_vertices [r][c][122] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][123] = 1.0f;
                filaments_vertices [r][c][124] = 0.0f;
                filaments_vertices [r][c][125] = 0.0f;
                filaments_vertices [r][c][126] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][127] = -1.0f;
                filaments_vertices [r][c][128] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][129] = 1.0f;
                filaments_vertices [r][c][130] = 0.0f;
                filaments_vertices [r][c][131] = 0.0f;
                filaments_vertices [r][c][132] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][133] = -1.0f;
                filaments_vertices [r][c][134] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][135] = 1.0f;
                filaments_vertices [r][c][136] = 0.0f;
                filaments_vertices [r][c][137] = 0.0f;
                filaments_vertices [r][c][138] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][139] = -1.0f;
                filaments_vertices [r][c][140] = 1.0f - delta_prime - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][141] = 1.0f;
                filaments_vertices [r][c][142] = 0.0f;
                filaments_vertices [r][c][143] = 0.0f;
                
                filaments_vertices [r][c][144] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime; // right
                filaments_vertices [r][c][145] = -1.0f;
                filaments_vertices [r][c][146] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][147] = 1.0f;
                filaments_vertices [r][c][148] = 0.0f;
                filaments_vertices [r][c][149] = 0.0f;
                filaments_vertices [r][c][150] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][151] = -1.0f;
                filaments_vertices [r][c][152] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][153] = 1.0f;
                filaments_vertices [r][c][154] = 0.0f;
                filaments_vertices [r][c][155] = 0.0f;
                filaments_vertices [r][c][156] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][157] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][158] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][159] = 1.0f;
                filaments_vertices [r][c][160] = 0.0f;
                filaments_vertices [r][c][161] = 0.0f;
                filaments_vertices [r][c][162] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][163] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][164] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][165] = 1.0f;
                filaments_vertices [r][c][166] = 0.0f;
                filaments_vertices [r][c][167] = 0.0f;
                filaments_vertices [r][c][168] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][169] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][170] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][171] = 1.0f;
                filaments_vertices [r][c][172] = 0.0f;
                filaments_vertices [r][c][173] = 0.0f;
                filaments_vertices [r][c][174] = -0.5*L + (float(c)*(delta_prime+R)) + delta_prime;
                filaments_vertices [r][c][175] = -1.0f;
                filaments_vertices [r][c][176] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][177] = 1.0f;
                filaments_vertices [r][c][178] = 0.0f;
                filaments_vertices [r][c][179] = 0.0f;
                
                filaments_vertices [r][c][180] = -0.5*L + (float(c)*(delta_prime+R)); // left
                filaments_vertices [r][c][181] = -1.0f;
                filaments_vertices [r][c][182] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][183] = 1.0f;
                filaments_vertices [r][c][184] = 0.0f;
                filaments_vertices [r][c][185] = 0.0f;
                filaments_vertices [r][c][186] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][187] = -1.0f;
                filaments_vertices [r][c][188] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][189] = 1.0f;
                filaments_vertices [r][c][190] = 0.0f;
                filaments_vertices [r][c][191] = 0.0f;
                filaments_vertices [r][c][192] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][193] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][194] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][195] = 1.0f;
                filaments_vertices [r][c][196] = 0.0f;
                filaments_vertices [r][c][197] = 0.0f;
                filaments_vertices [r][c][198] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][199] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][200] = 1.0f - (float(r)*(delta_prime+R)) - delta_prime;
                filaments_vertices [r][c][201] = 1.0f;
                filaments_vertices [r][c][202] = 0.0f;
                filaments_vertices [r][c][203] = 0.0f;
                filaments_vertices [r][c][204] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][205] = -1.0f + Num_subunits_array[r][c][i]*delta_prime + (float(r)/float(m))*delta_prime + dz_bend_array[r][c][i]/opengl_window_scaling_metric + dz_elas_array[r][c][i]/opengl_window_scaling_metric;
                filaments_vertices [r][c][206] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][207] = 1.0f;
                filaments_vertices [r][c][208] = 0.0f;
                filaments_vertices [r][c][209] = 0.0f;
                filaments_vertices [r][c][210] = -0.5*L + (float(c)*(delta_prime+R));
                filaments_vertices [r][c][211] = -1.0f;
                filaments_vertices [r][c][212] = 1.0f - (float(r)*(delta_prime+R));
                filaments_vertices [r][c][213] = 1.0f;
                filaments_vertices [r][c][214] = 0.0f;
                filaments_vertices [r][c][215] = 0.0f;
                
            }
            
            
        }}
            

//        GLfloat filaments_vertices [n_filaments][24];
//        GLuint filaments_indices [n_filaments][6];
//            
//            
//        for (int k=0; k<n_filaments; k++){
//                            
//            if ( k < 0.5*(n_filaments - n_central) || k >= 0.5*(n_filaments + n_central) ){
//                    
//                filaments_vertices [k][0] = -0.5*L + (float(k)*(delta_prime+R)) + delta_prime; // Top right x
//                filaments_vertices [k][1] = -0.75f + Num_subunits_array[k][i]*delta_prime + (float(k)/float(n_filaments))*delta_prime + dz_bend_array[k][i]/opengl_window_scaling_metric + dz_elas_array[k][i]/opengl_window_scaling_metric; // Top right y
//                filaments_vertices [k][2] = 0.0f;
//                filaments_vertices [k][3] = 1.0f;
//                filaments_vertices [k][4] = 0.0f;
//                filaments_vertices [k][5] = 0.0f;
//                filaments_vertices [k][6] = -0.5*L + (float(k)*(delta_prime+R)) + delta_prime; // Bottom right x
//                filaments_vertices [k][7] = -0.75f; // Bottom right y
//                filaments_vertices [k][8] = 0.0f;
//                filaments_vertices [k][9] = 1.0f;
//                filaments_vertices [k][10] = 0.0f;
//                filaments_vertices [k][11] = 0.0f;
//                filaments_vertices [k][12] = -0.5*L + (float(k)*(delta_prime+R)); // Bottom Left x
//                filaments_vertices [k][13] = -0.75f; // Bottom left y
//                filaments_vertices [k][14] = 0.0f;
//                filaments_vertices [k][15] = 1.0f;
//                filaments_vertices [k][16] = 0.0f;
//                filaments_vertices [k][17] = 0.0f;
//                filaments_vertices [k][18] = -0.5*L + (float(k)*(delta_prime+R)); // Top left x
//                filaments_vertices [k][19] = -0.75f + Num_subunits_array[k][i]*delta_prime + (float(k)/float(n_filaments))*delta_prime + dz_bend_array[k][i]/opengl_window_scaling_metric + dz_elas_array[k][i]/opengl_window_scaling_metric; // Top left y
//                filaments_vertices [k][20] = 0.0f;
//                filaments_vertices [k][21] = 1.0f;
//                filaments_vertices [k][22] = 0.0f;
//                filaments_vertices [k][23] = 0.0f;
//                }
//                else {
//                    filaments_vertices [k][0] = -0.5*L + (float(k)*(delta_prime+R)) + delta_prime; // Top right x
//                    filaments_vertices [k][1] = -0.75f + Num_subunits_array[k][i]*delta_prime + (float(k)/float(n_filaments))*delta_prime + dz_bend_array[k][i]/opengl_window_scaling_metric + dz_elas_array[k][i]/opengl_window_scaling_metric; // Top right y
//                    filaments_vertices [k][2] = 0.0f;
//                    filaments_vertices [k][3] = 0.0f;
//                    filaments_vertices [k][4] = 1.0f;
//                    filaments_vertices [k][5] = 0.0f;
//                    filaments_vertices [k][6] = -0.5*L + (float(k)*(delta_prime+R)) + delta_prime; // Bottom right x
//                    filaments_vertices [k][7] = -0.75f; // Bottom right y
//                    filaments_vertices [k][8] = 0.0f;
//                    filaments_vertices [k][9] = 0.0f;
//                    filaments_vertices [k][10] = 1.0f;
//                    filaments_vertices [k][11] = 0.0f;
//                    filaments_vertices [k][12] = -0.5*L + (float(k)*(delta_prime+R)); // Bottom Left x
//                    filaments_vertices [k][13] = -0.75f; // Bottom left y
//                    filaments_vertices [k][14] = 0.0f;
//                    filaments_vertices [k][15] = 0.0f;
//                    filaments_vertices [k][16] = 1.0f;
//                    filaments_vertices [k][17] = 0.0f;
//                    filaments_vertices [k][18] = -0.5*L + (float(k)*(delta_prime+R)); // Top left x
//                    filaments_vertices [k][19] = -0.75f + Num_subunits_array[k][i]*delta_prime + (float(k)/float(n_filaments))*delta_prime + dz_bend_array[k][i]/opengl_window_scaling_metric + dz_elas_array[k][i]/opengl_window_scaling_metric; // Top left y
//                    filaments_vertices [k][20] = 0.0f;
//                    filaments_vertices [k][21] = 0.0f;
//                    filaments_vertices [k][22] = 1.0f;
//                    filaments_vertices [k][23] = 0.0f;
//                    }
//                }
//            
//            
//        for (int k=0; k<n_filaments; k++){
//            filaments_indices [k][0] = 0;
//            filaments_indices [k][1] = 1;
//            filaments_indices [k][2] = 2;
//            filaments_indices [k][3] = 2;
//            filaments_indices [k][4] = 3;
//            filaments_indices [k][5] = 0;}
        

            
        GLfloat test_filaments [216];
            
        for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
            
            for (int j=0; j<216; j++){test_filaments [j] = filaments_vertices[r][c][j];}

            glBindVertexArray(VAOs_for_the_filaments[c + m*r]);
            
            glBindBuffer(GL_ARRAY_BUFFER, VBOs_for_the_filaments[c + m*r]);
            glBufferData(GL_ARRAY_BUFFER, sizeof(test_filaments), test_filaments, GL_STREAM_DRAW);

            glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * ) 0 );
            glEnableVertexAttribArray( 0 );
            glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
            glEnableVertexAttribArray( 1 );
            
            glBindVertexArray( 0 );
        }}
            
        
          
//
//                   GLfloat vertices[] =
//                 {
//                 -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, // back
//                 0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
//                 0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
//                 0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
//                 -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
//                 -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
//                 
//                 -0.5f, -0.5f,  0.5f,  0.0f, 0.0f, // front
//                 0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
//                 0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
//                 0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
//                 -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
//                 -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
//                 
//                 -0.5f,  0.5f,  0.5f,  1.0f, 0.0f, // left
//                 -0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
//                 -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
//                 -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
//                 -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
//                 -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
//                 
//                 0.5f,  0.5f,  0.5f,  1.0f, 0.0f, // right
//                 0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
//                 0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
//                 0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
//                 0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
//                 0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
//                 
//                 -0.5f, -0.5f, -0.5f,  0.0f, 1.0f, // bottom
//                 0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
//                 0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
//                 0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
//                 -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
//                 -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
//                 
//                 -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, // top
//                 0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
//                 0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
//                 0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
//                 -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
//                 -0.5f,  0.5f, -0.5f,  0.0f, 1.0f
//                 };
            
                
                
                //            GLfloat vertices[] =
                //            {
                //                -d, -d, -d,  0.0f, 0.0f, // back
                //                d, -d, -d,  1.0f, 0.0f,
                //                d, d, -d,  1.0f, 1.0f,
                //                d, d, -d,  1.0f, 1.0f,
                //                -d, d, -d,  0.0f, 1.0f,
                //                -d, -d, -d,  0.0f, 0.0f,
                //
                //                -d, -d, d,  0.0f, 0.0f, // back
                //                d, -d, d,  1.0f, 0.0f,
                //                d, d, d,  1.0f, 1.0f,
                //                d, d, d,  1.0f, 1.0f,
                //                -d, d, d,  0.0f, 1.0f,
                //                -d, -d, d,  0.0f, 0.0f,
                //
                //                -d, d, d,  1.0f, 0.0f, // left
                //                -d, d, -d,  1.0f, 1.0f,
                //                -d, -d, -d,  0.0f, 1.0f,
                //                -d, -d, -d,  0.0f, 1.0f,
                //                -d, -d, d,  0.0f, 0.0f,
                //                -d, d, d,  1.0f, 0.0f,
                //
                //                d, d, d,  1.0f, 0.0f, // left
                //                d, d, -d,  1.0f, 1.0f,
                //                d, -d, -d,  0.0f, 1.0f,
                //                d, -d, -d,  0.0f, 1.0f,
                //                d, -d, d,  0.0f, 0.0f,
                //                d, d, d,  1.0f, 0.0f,
                //
                //                -d, -d, -d,  0.0f, 1.0f, // bottom
                //                d, -d, -d,  1.0f, 1.0f,
                //                d, -d, d,  1.0f, 0.0f,
                //                d, -d, d,  1.0f, 0.0f,
                //                -d, -d, d,  0.0f, 0.0f,
                //                -d, -d, -d,  0.0f, 1.0f,
                //
                //                -d, d, -d,  0.0f, 1.0f, // bottom
                //                d, d, -d,  1.0f, 1.0f,
                //                d, d, d,  1.0f, 0.0f,
                //                d, d, d,  1.0f, 0.0f,
                //                -d, d, d,  0.0f, 0.0f,
                //                -d, d, -d,  0.0f, 1.0f,
                //            };
                
                
                
                        glm::vec3 cubePositions[] =
                {
                    
                                      glm::vec3( 0.0f, 0.0f, 0.0f ),
                    //                glm::vec3( -3.8f, -2.0f, -12.3f ),
                    //                glm::vec3( 2.4f, -0.4f, -3.5f ),
                    //                glm::vec3( -1.7f, 3.0f, -7.5f ),
                    //                glm::vec3( 1.3f, -2.0f, -2.5f ),
                    //                glm::vec3( 1.5f, 2.0f, -2.5f ),
                    //                glm::vec3( 1.5f, 0.2f, -1.5f ),
                    //                glm::vec3( -1.3f, 1.0f, -1.5f )
                };
                
                
                // Bind our Vertex Array Object first, then bind and set our buffers and pointers.
                //            glBindVertexArray( VAO );
                //            
                //            glBindBuffer( GL_ARRAY_BUFFER, VBO );
                //            glBufferData( GL_ARRAY_BUFFER, sizeof( vertices ), vertices, GL_STATIC_DRAW );
                //            
                //            // Position attribute
                //            glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof( GLfloat ), ( GLvoid * )0 );
                //            glEnableVertexAttribArray( 0 );
                //            // TexCoord attribute
                //            glVertexAttribPointer( 2, 2, GL_FLOAT, GL_FALSE, 5 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
                //            glEnableVertexAttribArray( 2 );
                //            
                //            glBindVertexArray( 0 ); // Unbind VAO
            
                
                
            
        ourShader.Use();
            
            
            for (int j=0; j<=5000; j++){
                
                GLfloat point[] = { 0.0f, 0.0f, 0.0f,  0.5f, 0.5f, 0.0f };
                
                glBindVertexArray(VAO_for_a_point);
                glBindBuffer(GL_ARRAY_BUFFER, VBO_for_a_point);
                glBufferData(GL_ARRAY_BUFFER, sizeof(point), point, GL_STATIC_DRAW);
                glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * ) 0 );
                glEnableVertexAttribArray( 0 );
                glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
                glEnableVertexAttribArray( 1 );
                glBindVertexArray( 0 );
                
                glBindVertexArray(VAO_for_a_point);
                glDrawArrays(GL_POINT, 0, 3);
                glBindVertexArray(0);
            }
            
                

            glBindVertexArray(VAO_for_the_mem);
            glDrawArrays( GL_TRIANGLES, 0, 36 );
            glBindVertexArray(0);


            
//        glBindVertexArray(VAO_for_the_mem[0]);
//        glDrawElements( GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0 ); // The last argument specifies the starting index of the vertex array we'd like to draw; we just leave this at 0. The second argument specifies how many vertices we want to draw.
//        glBindVertexArray(0);

            
            
        for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
            
            glBindVertexArray(VAOs_for_the_filaments[c + m*r]);
            glDrawArrays( GL_TRIANGLES, 0, 36 );
            glBindVertexArray(0);
            
        }}
            

            
            glm::mat4 projection;
            projection = glm::perspective(camera.GetZoom( ), (GLfloat)SCREEN_WIDTH/(GLfloat)SCREEN_HEIGHT, 0.1f, 1000.0f);
            
            // Create camera transformation
            glm::mat4 view;
            view = camera.GetViewMatrix( );
            
            // Get the uniform locations
            GLint modelLoc = glGetUniformLocation( ourShader.Program, "model" );
            GLint viewLoc = glGetUniformLocation( ourShader.Program, "view" );
            GLint projLoc = glGetUniformLocation( ourShader.Program, "projection" );
//
            // Pass the matrices to the shader
            glUniformMatrix4fv( viewLoc, 1, GL_FALSE, glm::value_ptr( view ) );
            glUniformMatrix4fv( projLoc, 1, GL_FALSE, glm::value_ptr( projection ) );
            
            
            glBindVertexArray( VAO );
            for( GLuint i = 0; i < 1; i++ )
            {
                // Calculate the model matrix for each object and pass it to shader before drawing
                glm::mat4 model;
                model = glm::translate( model, cubePositions[i] );
                GLfloat angle = 50.0f * i;
                model = glm::rotate(model, angle, glm::vec3( 1.0f, 0.3f, 0.5f ) );
                glUniformMatrix4fv( modelLoc, 1, GL_FALSE, glm::value_ptr( model ) );
                
               // glDrawArrays( GL_TRIANGLES, 0, 36 );
            }
            glBindVertexArray( 0 );
            
            
            
            glfwSwapBuffers( window ); // Swap the screen buffers
        }
    }
    
    glDeleteVertexArrays(1, &VAO_for_a_point);
    glDeleteVertexArrays(1, &VAO_for_the_mem);
    glDeleteVertexArrays(n_filaments, VAOs_for_the_filaments);

    glDeleteBuffers(1, &VBO_for_a_point);
    glDeleteBuffers(1, &VBO_for_the_mem);
    glDeleteBuffers(n_filaments, VBOs_for_the_filaments);

    glDeleteVertexArrays( 1, &VAO );
    glDeleteBuffers( 1, &VBO );
    
    glfwTerminate( );
    
    return EXIT_SUCCESS;
        
}


// Moves/alters the camera positions based on user input
void DoMovement( )
{
    // Camera controls
    if( keys[GLFW_KEY_W] || keys[GLFW_KEY_UP] )
    {
        camera.ProcessKeyboard( FORWARD, deltaTime );
    }
    
    if( keys[GLFW_KEY_S] || keys[GLFW_KEY_DOWN] )
    {
        camera.ProcessKeyboard( BACKWARD, deltaTime );
    }
    
    if( keys[GLFW_KEY_A] || keys[GLFW_KEY_LEFT] )
    {
        camera.ProcessKeyboard( LEFT, deltaTime );
    }
    
    if( keys[GLFW_KEY_D] || keys[GLFW_KEY_RIGHT] )
    {
        camera.ProcessKeyboard( RIGHT, deltaTime );
    }
}

// Is called whenever a key is pressed/released via GLFW
void KeyCallback( GLFWwindow *window, int key, int scancode, int action, int mode )
{
    if( key == GLFW_KEY_ESCAPE && action == GLFW_PRESS )
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    
    if ( key >= 0 && key < 1024 )
    {
        if( action == GLFW_PRESS )
        {
            keys[key] = true;
        }
        else if( action == GLFW_RELEASE )
        {
            keys[key] = false;
        }
    }
}

void MouseCallback( GLFWwindow *window, double xPos, double yPos )
{
    if( firstMouse )
    {
        lastX = xPos;
        lastY = yPos;
        firstMouse = false;
    }
    
    GLfloat xOffset = xPos - lastX;
    GLfloat yOffset = lastY - yPos;  // Reversed since y-coordinates go from bottom to left
    
    lastX = xPos;
    lastY = yPos;
    
    camera.ProcessMouseMovement( xOffset, yOffset );
}


void ScrollCallback( GLFWwindow *window, double xOffset, double yOffset )
{
    camera.ProcessMouseScroll( yOffset );
}

