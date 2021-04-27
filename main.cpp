#define _USE_MATH_DEFINES
#include <math.h>

#include "engine_reader.h"


using namespace std;

float px = 00, py = 30, pz = 40;
float radius = 400;
float lx = 0.0, ly = 0.0, lz = 0.0;
float alpha = 45.0, beta = 45.0;
int frame = 0, timebase = 0;
int startX, startY, tracking = 0;

float moveZ = 0.0;

struct Scene scene;
GLuint* buffers;
int *n_verteces;
int buffer_counter, draw_counter;
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   LETS FIND THE F*CKING ERROR                                      //
//                                                                                                    //

GLenum error;

static void GlClearError() {
    while (glGetError() != GL_NO_ERROR);
}

static bool GLLogCall() {

    while (error != glGetError()) {
        std::cout << "OpenGL Error: " << error << std::endl;
        return false;
    }
    return true;
}

//                                                                                                  //
//                                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////////


void draw_vbo() {
    glBindBuffer(GL_ARRAY_BUFFER, buffers[draw_counter]);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glDrawArrays(GL_TRIANGLES, 0, buffers[draw_counter++]);
}

void fill_buffers(vector<struct Group> groups) {

    for (int k = 0; k < groups.size(); k++) {

        fill_buffers(groups[k].child);

        vector<vector<struct Point>> models = groups[k].models;

        for (int i = 0; i < models.size(); i++) {

            n_verteces[buffer_counter] += models[i].size();

            glBindBuffer(GL_ARRAY_BUFFER, buffers[buffer_counter++]);

            vector<struct Point> model = models[i];

            float* verteces = (float*)malloc(sizeof(float*) *model.size() * 3);

            for (int j = 0; j < model.size(); j++) {

                verteces[3 * j] = model[j].x;
                verteces[3 * j + 1] = model[j].y;
                verteces[3 * j + 2] = model[j].z;
            }
            //printf("x = %f, y = %f, z = %f\n", verteces[3 * j], verteces[3 * j + 1], verteces[3 * j + 2]);

            glBufferData(GL_ARRAY_BUFFER, sizeof(float) * model.size() * 3, verteces, GL_STATIC_DRAW);

            free(verteces);
        }
    }
}

void prepare_vbo_data() {

    n_verteces = (int*)malloc(sizeof(int) * scene.nModels);

    glEnableClientState(GL_VERTEX_ARRAY);

    buffer_counter = 0;

    buffers = (GLuint*)malloc(sizeof(GLuint) * scene.nModels);

    glGenBuffers(scene.nModels, buffers);

    fill_buffers(scene.groups);
}


void refreshCam() {

    px = radius * cos(beta) * sin(alpha);
    py = radius * sin(beta);
    pz = radius * cos(beta) * cos(alpha);
}

void changeSize(int w, int h) {
    // Prevent a divide by zero, when window is too short
    // (you cant make a window with zero width).
    if (h == 0)
        h = 1;

    // compute window's aspect ratio
    double ratio = w * 1.0 / h;

    // Set the projection matrix as current
    glMatrixMode(GL_PROJECTION);
    // Load Identity Matrix
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);

    // Set perspective
    gluPerspective(45.0f, ratio, 1.0f, 1000.0f);

    // return to the model view matrix mode
    glMatrixMode(GL_MODELVIEW);
}


void draw_scene(vector<struct Group> groups, int time) {

    for (int i = 0; i < groups.size(); i++) {
        int k = 0;
        struct Group group = groups[i];
        glPushMatrix();
        {
            if (!group.colors.empty()) {
                draw_color(group.colors[k++]);
            }
            else {
                glColor3f(1.0, 1.0, 1.0);
            }
            draw_gt(group, time);
            draw_vbo();
            //draw_models(group);
            draw_scene(group.child, time);

        }
        glPopMatrix();
    }
}

void show_fps(int time) {
    float fps;
    char s[64];

    frame++;

    if (time - timebase > 1000) {
        fps = frame * 1000.0 / (time - timebase);
        timebase = time;
        frame = 0;
        sprintf(s, "FPS: %6.2f", fps);
        glutSetWindowTitle(s);
    }
}


void renderScene(void) {

    int time;
    draw_counter = 0;
        
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set the camera
    glLoadIdentity();
    gluLookAt(px, py, pz,
              lx, ly, lz + moveZ,
              0.0f, 1.0f, 0.0f);

    // put the geometric transformations here
    time = glutGet(GLUT_ELAPSED_TIME);

    draw_scene(scene.groups, time);

    time = glutGet(GLUT_ELAPSED_TIME);
    show_fps(time);


    // End of frame
    glutPostRedisplay();
    glutSwapBuffers();
}


void special_keyboard(int key_code, int a, int b) {
    switch (key_code) {

    case GLUT_KEY_UP:
        beta += 0.1f;
        if (beta > 1.5f)
            beta = 1.5f;
        break;
    case GLUT_KEY_DOWN:
        beta -= 0.1f;
        if (beta < -1.5f)
            beta = -1.5f;
        break;
    case GLUT_KEY_RIGHT:
        alpha -= 0.1;
        break;
    case GLUT_KEY_LEFT:
        alpha += 0.1;
        break;
    default: break;
    }

    refreshCam();
    glutPostRedisplay();
}

void keyboard(unsigned char key_code, int a, int b) {
    switch (key_code) {
    case 27:	exit(0); break;
    case 'w':	radius -= 2; break;
    case 's':	radius += 2; break;
    case 'a':   moveZ -= 2; break;
    case 'd':   moveZ += 2; break;
    case 'p':	glPolygonMode(GL_FRONT, GL_POINT); break;
    case 'l':	glPolygonMode(GL_FRONT, GL_LINE); break;
    case 'o':	glPolygonMode(GL_FRONT, GL_FILL); break;

    default: break;
    }

    refreshCam();
    glutPostRedisplay();
}

void processMouseButtons(int button, int state, int xx, int yy) {

    if (state == GLUT_DOWN) {
        startX = xx;
        startY = yy;
        if (button == GLUT_LEFT_BUTTON)
            tracking = 1;
        else if (button == GLUT_RIGHT_BUTTON)
            tracking = 2;
        else
            tracking = 0;
    }
    else if (state == GLUT_UP) {
        if (tracking == 1) {
            alpha += (xx - startX);
            beta += (yy - startY);
        }
        else if (tracking == 2) {

            radius -= yy - startY;
            if (radius < 3)
                radius = 3.0;
        }
        tracking = 0;
    }

}

void processMouseMotion(int xx, int yy) {

    int deltaX, deltaY;
    int alphaAux, betaAux;
    int rAux;

    if (!tracking)
        return;

    deltaX = xx - startX;
    deltaY = yy - startY;

    if (tracking == 1) {


        alphaAux = alpha + deltaX;
        betaAux = beta + deltaY;

        if (betaAux > 85.0)
            betaAux = 85.0;
        else if (betaAux < -85.0)
            betaAux = -85.0;

        rAux = radius;
    }
    else if (tracking == 2) {

        alphaAux = alpha;
        betaAux = beta;
        rAux = radius - deltaY;
        if (rAux < 3)
            rAux = 3;
    }
    px = rAux * sin(alphaAux * M_PI / 180.0) * cos(betaAux * M_PI / 180.0);
    pz = rAux * cos(alphaAux * M_PI / 180.0) * cos(betaAux * M_PI / 180.0);
    py = rAux * sin(betaAux * M_PI / 180.0);
}




int main(int argc, char **argv){

    if (argc < 2) {
        printf("\nYUP, you have an Input Error   [main.cpp -> (Line %d) ]\n\n", __LINE__);
        return -1;
    }

    string fileDir = "../../Models/";
    string xmlFile = fileDir + argv[1];
    TiXmlDocument doc(xmlFile.c_str());
    if (!doc.LoadFile()){

        printf("ERROR: Error opening .xml File!!\n");
        return -2;
    }
    TiXmlElement *root = doc.RootElement();
    for (TiXmlElement *group = root->FirstChild("group")->ToElement(); group; group = group->NextSiblingElement()) {
        load_scene(&scene, group);
    }

    // init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(1000, 800);
    glutCreateWindow("Practical Assignment CG - Phase 3");

    // Required callback registry
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);

    // put here the registration of the keyboard callbacks
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special_keyboard);
    glutMouseFunc(processMouseButtons);
    glutMotionFunc(processMouseMotion);

    //  OpenGL settings
    glPolygonMode(GL_FRONT, GL_LINE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    if (glewInit() != GLEW_OK) {
        printf("Error initializating glew\n");
        return -1;
    }

    refreshCam();

    prepare_vbo_data();

    // enter GLUT's main cycle
    glutMainLoop();

    return 1;
}





