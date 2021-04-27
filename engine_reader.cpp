#include "engine_reader.h"

struct Point init_point(float x, float y, float z)
{
    struct Point p;
    p.x = x;
    p.y = y;
    p.z = z;
    return p;
}

struct GT init_translate()
{
    struct GT t;
    t.p = init_point(0, 0, 0);
    t.gt_type = TRANSLATE;

    return t;
} 

struct GT init_scale(void)
{
    struct GT s;
    s.p = init_point(1, 1, 1);
    s.gt_type = SCALE;

    return s;
}

struct GT init_rotate(void) {
    struct GT r;
    r.p = init_point(0, 0, 0);
    r.gt_type = ROTATE;
    r.r_angle = 0.0;
    r.rotate_time = 0.0;

    return r;
}

struct GT init_translate_anim(void) {
    struct GT ta;
    ta.p = init_point(0, 0, 0);
    ta.gt_type = TRANSLATE_ANIM;
    ta.translate_time = 0;

    return ta;
}

struct GT init_rotate_anim(void)
{
    struct GT ra;
    ra.p = init_point(0, 0, 0);
    ra.gt_type = ROTATE_ANIM;
    ra.rotate_time = 0;

    return ra;
}

void buildRotMatrix(float* x, float* y, float* z, float* m) {
    m[0] = x[0];    m[4] = y[0];    m[8] = z[0];    m[12] = 0;
    m[1] = x[1];    m[5] = y[1];    m[9] = z[1];    m[13] = 0;
    m[2] = x[2];    m[6] = y[2];    m[10] = z[2];   m[14] = 0;
    m[3] = 0;       m[7] = 0;       m[11] = 0;      m[15] = 1;

}

void cross(float* a, float* b, float* res){

    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}
  
void normalize(float* a) {

    float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] = a[0] / l;
    a[1] = a[1] / l;
    a[2] = a[2] / l;
}

float length(float* v) {

    float res = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return res;
}

void multMatrixVector(float* m, float* v, float* res) {

    for (int j = 0; j < 4; ++j) {
        res[j] = 0;
        for (int k = 0; k < 4; ++k) {
            res[j] += v[k] * m[j * 4 + k];
        }
    }
}

//Function that loads the translation of the .xml file
void load_translate(TiXmlElement* transf, struct Group* g)
{
    struct GT point = init_translate();

    if (transf->Attribute("X"))
        point.p.x = atof(transf->Attribute("X"));
    if (transf->Attribute("Y"))
        point.p.y = atof(transf->Attribute("Y"));
    if (transf->Attribute("Z"))
        point.p.z = atof(transf->Attribute("Z"));

    g->gt.push_back(point);
}

//Function that loads the scale data from the .xml file
void load_scale(TiXmlElement* transf, struct Group* g) {
    struct GT scale = init_scale();

    if (transf->Attribute("X"))
        scale.p.x = atof(transf->Attribute("X"));
    if (transf->Attribute("Y"))
        scale.p.y = atof(transf->Attribute("Y"));
    if (transf->Attribute("Z"))
        scale.p.z = atof(transf->Attribute("Z"));

    g->gt.push_back(scale);
}

//Function that loads the rotation data from the .xml file
void load_rotate(TiXmlElement* transf, struct Group* g) {
    struct GT rotate = init_rotate();
    
    if (transf->Attribute("time"))
        rotate.rotate_time = atof(transf->Attribute("time"));
    if (transf->Attribute("angle"))
        rotate.r_angle = atof(transf->Attribute("angle"));
    if (transf->Attribute("axisX"))
        rotate.p.x = atof(transf->Attribute("axisX"));
    if (transf->Attribute("axisY"))
        rotate.p.y = atof(transf->Attribute("axisY"));
    if (transf->Attribute("axisZ"))
        rotate.p.z = atof(transf->Attribute("axisZ"));

    g->gt.push_back(rotate);
}

void load_rotate_anim(TiXmlElement* transf, struct Group* g)
{
    struct GT ra = init_rotate_anim();

    ra.rotate_time = atof(transf->Attribute("time")) * 100;

    if (transf->Attribute("axisX"))
        ra.p.x = atof(transf->Attribute("axisX"));
    if (transf->Attribute("axisY"))
        ra.p.y = atof(transf->Attribute("axisY"));
    if (transf->Attribute("axisZ"))
        ra.p.z = atof(transf->Attribute("axisZ"));

    g->gt.push_back(ra);
}

void load_translate_anim(TiXmlElement* transf, struct Group* g)
{
    struct GT ta = init_translate_anim();

    ta.translate_time = atof(transf->Attribute("time")) * 1000;

    if (ta.translate_time > 0) {

        for (TiXmlElement* points = transf->FirstChildElement("point"); points; points = points->NextSiblingElement()) {

            Point pp = init_point(0, 0, 0);
            if (points->Attribute("X"))
                pp.x = atof(points->Attribute("X"));
            if (points->Attribute("Y"))
                pp.y = atof(points->Attribute("Y"));
            if (points->Attribute("Z"))
                pp.z = atof(points->Attribute("Z"));

            ta.control_points.push_back(pp);
        }
    }
    g->gt.push_back(ta);
}

void load_color(TiXmlElement* color, struct Group* g) {

    struct Point col = init_point(0, 0, 0);

    if (color->Attribute("R"))
        col.x = atof(color->Attribute("R"));
    if (color->Attribute("G"))
        col.y = atof(color->Attribute("G"));
    if (color->Attribute("B"))
        col.z = atof(color->Attribute("B"));

    g->colors.push_back(col);
}

//Function that loads the models data from the .xml and subsequently the .3d file
void load_models(TiXmlElement* models, struct Group* g) {

    if (TiXmlElement* color = models->FirstChildElement("color")) {
        load_color(color, g);
    }


    for (TiXmlElement* model = models->FirstChildElement("model"); model; model = model->NextSiblingElement()) {
        const char* file = model->Attribute("file");

        if (file) {

            string fileDir = "../../Models/";
            fstream f;
            string file3d = fileDir + file;
            f.open(file3d, ios::in);

            if (f.is_open()) {

                string line;
                vector<struct Point> points;
                getline(f, line);

                while (!line.empty()) {

                    double x, y, z;

                    sscanf(line.c_str(), "%lf %lf %lf\n", &x, &y, &z);

                    struct Point p = init_point(x, y, z);
                    points.push_back(p);
                    getline(f, line);
                }

                g->models.push_back(points);
                f.close();
            }
            else {
                printf("\nYUP, you have an error opening 3D file (Line %d) \n\n", __LINE__);
            }
        }
        else {
            printf("\nYUP, you have an error opening XML file (Line %d) \n\n", __LINE__);
        }
    }
}

struct Group process_groups(TiXmlElement* group, struct Scene* scene)
{
    struct Group g;

    if (!group) {
        return g;
    }

    for (TiXmlElement* childGroup = group->FirstChildElement("group"); childGroup; childGroup = childGroup->NextSiblingElement()) {

        g.child.push_back(process_groups(childGroup, scene));
    }

    TiXmlElement* trans = group->FirstChildElement("translate");
    if (trans) {

        if (trans->Attribute("time")) {

            load_translate_anim(trans, &g);
        }
        else {
            load_translate(trans, &g);
        }
    }

    TiXmlElement* scale = group->FirstChildElement("scale");
    if (scale) {
        load_scale(scale, &g);
    }

    TiXmlElement* rotate = group->FirstChildElement("rotate");
    if (rotate) {

        if (rotate->Attribute("time")){

            load_rotate_anim(rotate, &g);
        }
        else{

            load_rotate(rotate, &g);
        }
    }

    TiXmlElement* models = group->FirstChildElement("models");
    if (models){

        scene->nModels++;

        load_models(models, &g);
    }
    return g;
}

void load_scene(struct Scene* scene, TiXmlElement* group)
{
    scene->groups.push_back(process_groups(group, scene));
}

void draw_models(struct Group g)
{
    for (unsigned i = 0; i < g.models.size(); i++)
    {
        vector<struct Point> model = g.models[i];
        for (unsigned j = 0; j < model.size(); j += 3)
        {
            glBegin(GL_TRIANGLES);
            glVertex3f(model[j].x, model[j].y, model[j].z);
            glVertex3f(model[j + 1].x, model[j + 1].y, model[j + 1].z);
            glVertex3f(model[j + 2].x, model[j + 2].y, model[j + 2].z);
            glEnd();
        }
    }
}

void getCatmullRomPoint(float t, struct Point p0, struct Point p1, struct Point p2, struct Point p3, float* res, float* deriv)
{
    // catmull-rom matrix
    float m[16] = { -0.5f, 1.5f, -1.5f, 0.5f,
                   1.0f, -2.5f, 2.0f, -0.5f,
                   -0.5f, 0.0f, 0.5f, 0.0f,
                   0.0f, 1.0f, 0.0f, 0.0f };

    // reset res and deriv
    res[0] = 0.0;
    res[1] = 0.0;
    res[2] = 0.0;
    deriv[0] = 0.0;
    deriv[1] = 0.0;
    deriv[2] = 0.0;
    // Compute A = M * P
    float Ax[4], Ay[4], Az[4];

    float Px[4], Py[4], Pz[4];

    Px[0] = p0.x;
    Px[1] = p1.x;
    Px[2] = p2.x;
    Px[3] = p3.x;

    Py[0] = p0.y;
    Py[1] = p1.y;
    Py[2] = p2.y;
    Py[3] = p3.y;

    Pz[0] = p0.z;
    Pz[1] = p1.z;
    Pz[2] = p2.z;
    Pz[3] = p3.z;

    multMatrixVector(m, Px, Ax);
    multMatrixVector(m, Py, Ay);
    multMatrixVector(m, Pz, Az);

    // Compute point res = T *A
    //  res = T * A
    float P[4] = { t * t * t, t * t, t, 1 };

    res[0] = P[0] * Ax[0] + P[1] * Ax[1] + P[2] * Ax[2] + P[3] * Ax[3];
    res[1] = P[0] * Ay[0] + P[1] * Ay[1] + P[2] * Ay[2] + P[3] * Ay[3];
    res[2] = P[0] * Az[0] + P[1] * Az[1] + P[2] * Az[2] + P[3] * Az[3];

    // deriv = T' * A
    float Pl[4] = { 3 * (t * t), 2 * t, 1, 0 };

    deriv[0] = Pl[0] * Ax[0] + Pl[1] * Ax[1] + Pl[2] * Ax[2] + Pl[3] * Ax[3];
    deriv[1] = Pl[0] * Ay[0] + Pl[1] * Ay[1] + Pl[2] * Ay[2] + Pl[3] * Ay[3];
    deriv[2] = Pl[0] * Az[0] + Pl[1] * Az[1] + Pl[2] * Az[2] + Pl[3] * Az[3];
}

void getGlobalCatmullRomPoint(float gt, float* pos, float* deriv, vector<Point> controlPoints)
{
    int pointCount = controlPoints.size();
    float t = gt * pointCount; // this is the real global t
    int index = floor(t);      // which segment
    t = t - index;             // where within  the segment

    // indices store the points
    int indices[4];
    indices[0] = (index + pointCount - 1) % pointCount;
    indices[1] = (indices[0] + 1) % pointCount;
    indices[2] = (indices[1] + 1) % pointCount;
    indices[3] = (indices[2] + 1) % pointCount;

    getCatmullRomPoint(t, controlPoints[indices[0]], controlPoints[indices[1]], controlPoints[indices[2]], controlPoints[indices[3]], pos, deriv);
}

void renderCarmullRomCurve(struct GT gt)
{
    glBegin(GL_LINE_LOOP);
    float t = 0.0f;
    for (int i = 0; i < 1000; i++)
    {
        float P[3];
        float Pl[3];

        getGlobalCatmullRomPoint(t, (float*)P, (float*)Pl, gt.control_points);
        glVertex3f(P[0], P[1], P[2]);
        t += 0.001f;
    }
    glEnd();
}

void draw_gt(struct Group g, int elapsed)
{
    float Y[3] = { 0, 1, 0 };
    for (int i = 0; i < g.gt.size(); i++) {

        struct GT transf = g.gt[i];

        if (transf.gt_type == SCALE){ 
            glScalef(transf.p.x, transf.p.y, transf.p.z);
        }
        if (transf.gt_type == TRANSLATE){ 
            glTranslatef(transf.p.x, transf.p.y, transf.p.z);
        }
        if (transf.gt_type == ROTATE){ 
            glRotatef(transf.r_angle, transf.p.x, transf.p.y, transf.p.z);
        }
        if (transf.gt_type == ROTATE_ANIM){ 
            float angle = (360.0 * elapsed) / transf.rotate_time;
            glRotatef(angle, transf.p.x, transf.p.y, transf.p.z);
        }
        if (transf.gt_type == TRANSLATE_ANIM){
            float Y[3] = { 0, 1, 0 };
            float Z[3];
            float P[3];
            float Pl[3];
            float m[16];

            float t = (float)elapsed / (float)transf.translate_time;

            renderCarmullRomCurve(transf);

            getGlobalCatmullRomPoint(t, (float*)P, (float*)Pl, transf.control_points);
            normalize(Pl);

            // Z = P' x Y(i-1)
            cross(Pl, Y, Z);
            normalize(Z);

            // Y = Z x P'
            cross(Z, Pl, Y);
            normalize(Y);

            buildRotMatrix(Pl, Y, Z, m);

            glTranslatef(P[0], P[1], P[2]);

            glMultMatrixf(m);
        }
    }
}

void draw_color(struct Point color)         
{
    glColor3f(color.x, color.y, color.z);
}