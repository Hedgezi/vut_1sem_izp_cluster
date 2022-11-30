/**
 * Kostra programu pro 2. projekt IZP 2022/23
 *
 * Jednoducha shlukova analyza: 2D nejblizsi soused.
 * Single linkage
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sqrtf
#include <limits.h> // INT_MAX

#include <string.h> // memcpy
#include <time.h> // time

/*****************************************************************
 * Ladici makra. Vypnout jejich efekt lze definici makra
 * NDEBUG, napr.:
 *   a) pri prekladu argumentem prekladaci -DNDEBUG
 *   b) v souboru (na radek pred #include <assert.h>
 *      #define NDEBUG
 */
#ifdef NDEBUG
#define debug(s)
#define dfmt(s, ...)
#define dint(i)
#define dfloat(f)
#else

// vypise ladici retezec
#define debug(s) printf("- %s\n", s)

// vypise formatovany ladici vystup - pouziti podobne jako printf
#define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)

// vypise ladici informaci o promenne - pouziti dint(identifikator_promenne)
#define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)

// vypise ladici informaci o promenne typu float - pouziti
// dfloat(identifikator_promenne)
#define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)

#endif

#define checkmalloc(p) if ((p) == NULL) {fprintf(stderr, "can't allocate (on line %d).\n", __LINE__); return -1;}
#define checkpointer(p) if ((p) == NULL) {fprintf(stderr, "was given pointer to NULL (on line %d).\n", __LINE__); return -1;}
#define checkforerror(statement, error) if ((statement)) {fprintf(stderr, "error: " error); return -1;}
#define checkforerrorandfree(statement, error, n) if ((statement)) {fprintf(stderr, "error: " error); free_clusters(n); return -1;}


/*****************************************************************
 * Deklarace potrebnych datovych typu:
 *
 * TYTO DEKLARACE NEMENTE
 *
 *   struct obj_t - struktura objektu: identifikator a souradnice
 *   struct cluster_t - shluk objektu:
 *      pocet objektu ve shluku,
 *      kapacita shluku (pocet objektu, pro ktere je rezervovano
 *          misto v poli),
 *      ukazatel na pole shluku.
 */

struct obj_t {
    int id;
    float x;
    float y;
};

struct cluster_t {
    int size;
    int capacity;
    struct obj_t *obj;
};

/*****************************************************************
 * Deklarace potrebnych funkci.
 *
 * PROTOTYPY FUNKCI NEMENTE
 *
 * IMPLEMENTUJTE POUZE FUNKCE NA MISTECH OZNACENYCH 'TODO'
 *
 */

/*
 Inicializace shluku 'c'. Alokuje pamet pro cap objektu (kapacitu).
 Ukazatel NULL u pole objektu znamena kapacitu 0.
*/
int init_cluster(struct cluster_t *c, int cap)
{
    assert(c != NULL);
    assert(cap >= 0);

    c->size = 0;
    c->capacity = cap;
    if (cap > 0) {
        c->obj = malloc(cap*sizeof(struct obj_t));
        checkmalloc(c->obj);
    }
    else if (cap == 0)
        c->obj = NULL;
    else {
        fprintf(stderr, "error: negative capacity");
        exit(1);
        return -1;
    }
    return 0;
}

/*
 Odstraneni vsech objektu shluku a inicializace na prazdny shluk.
 */
int clear_cluster(struct cluster_t *c)
{
    assert(c != NULL);

    checkpointer(c);
    checkpointer(c->obj);

    free(c->obj);
    c->obj = NULL;
    c->capacity = 0;
    c->size = 0;
    return 0;
}

/// Chunk of cluster objects. Value recommended for reallocation.
const int CLUSTER_CHUNK = 10;

/*
 Zmena kapacity shluku 'c' na kapacitu 'new_cap'.
 */
struct cluster_t *resize_cluster(struct cluster_t *c, int new_cap)
{
    // TUTO FUNKCI NEMENTE
    assert(c);
    assert(c->capacity >= 0);
    assert(new_cap >= 0);

    if (c->capacity >= new_cap)
        return c;

    size_t size = sizeof(struct obj_t) * new_cap;

    void *arr = realloc(c->obj, size);
    if (arr == NULL)
        return NULL;

    c->obj = (struct obj_t*)arr;
    c->capacity = new_cap;
    return c;
}

/*
 Prida objekt 'obj' na konec shluku 'c'. Rozsiri shluk, pokud se do nej objekt
 nevejde.
 */
int append_cluster(struct cluster_t *c, struct obj_t obj)
{
    assert(c != NULL);

    checkpointer(c);

    if (c->size == c->capacity) {
        c = resize_cluster(c, c->capacity+CLUSTER_CHUNK);
        checkmalloc(c);
    }
    c->obj[c->size] = obj;
    c->size++;
    return 0;
}

/*
 Seradi objekty ve shluku 'c' vzestupne podle jejich identifikacniho cisla.
 */
void sort_cluster(struct cluster_t *c);

/*
 Do shluku 'c1' prida objekty 'c2'. Shluk 'c1' bude v pripade nutnosti rozsiren.
 Objekty ve shluku 'c1' budou serazeny vzestupne podle identifikacniho cisla.
 Shluk 'c2' bude nezmenen.
 */
int merge_clusters(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c2 != NULL);

    checkpointer(c1);
    checkpointer(c2);

    if (c1->capacity < c1->size+c2->size) {
        c1 = resize_cluster(c1, c1->size+c2->size);
    }
    for (int i = 0; i < c2->size; i++) {
        append_cluster(c1, c2->obj[i]);
    }
    return 0;
}

/**********************************************************************/
/* Prace s polem shluku */

/*
 Odstrani shluk z pole shluku 'carr'. Pole shluku obsahuje 'narr' polozek
 (shluku). Shluk pro odstraneni se nachazi na indexu 'idx'. Funkce vraci novy
 pocet shluku v poli.
*/
int remove_cluster(struct cluster_t *carr, int narr, int idx)
{
    assert(idx < narr);
    assert(narr > 0);

    free(carr[idx].obj);
    for (int i = idx; i < narr-1; i++) {
        carr[i] = carr[i+1];
    }
    return narr - 1;
}

/*
 Pocita Euklidovskou vzdalenost mezi dvema objekty.
 */
float obj_distance(struct obj_t *o1, struct obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    return sqrt(pow((o1->x - o2->x), 2)+pow((o1->y - o2->y), 2));
}

/*
 Pocita vzdalenost dvou shluku.
*/
float cluster_distance(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float dist = obj_distance(&c1->obj[0], &c2->obj[0]);

    for (int i = 0; i < c1->size; i++) {
        for (int j = 0; j < c2->size; j++) {
            dist = (dist > obj_distance(&c1->obj[i], &c2->obj[j])) ? obj_distance(&c1->obj[i], &c2->obj[j]) : dist;
        }
    }
    return dist;
}

/*
 Funkce najde dva nejblizsi shluky. V poli shluku 'carr' o velikosti 'narr'
 hleda dva nejblizsi shluky. Nalezene shluky identifikuje jejich indexy v poli
 'carr'. Funkce nalezene shluky (indexy do pole 'carr') uklada do pameti na
 adresu 'c1' resp. 'c2'.
*/
void find_neighbours(struct cluster_t *carr, int narr, int *c1, int *c2, float (*cd)(struct cluster_t *, struct cluster_t *))
{
    assert(narr > 0);

    float mindist = cd(&carr[0], &carr[1]);
    *c1 = 0; *c2 = 1;

    for (int i = 0; i < narr; i++) {
        for (int j = i+1; j < narr; j++) {
            if (mindist > cd(&carr[i], &carr[j])) {
                mindist = cd(&carr[i], &carr[j]);
                *c1 = i;
                *c2 = j;
            }
        }
    }
}

// pomocna funkce pro razeni shluku
static int obj_sort_compar(const void *a, const void *b)
{
    // TUTO FUNKCI NEMENTE
    const struct obj_t *o1 = (const struct obj_t *)a;
    const struct obj_t *o2 = (const struct obj_t *)b;
    if (o1->id < o2->id) return -1;
    if (o1->id > o2->id) return 1;
    return 0;
}

/*
 Razeni objektu ve shluku vzestupne podle jejich identifikatoru.
*/
void sort_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    qsort(c->obj, c->size, sizeof(struct obj_t), &obj_sort_compar);
}

/*
 Tisk shluku 'c' na stdout.
*/
void print_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    for (int i = 0; i < c->size; i++)
    {
        if (i) putchar(' ');
        printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
    }
    putchar('\n');
}

/*
Free all memory allocated for clusters and objects
*/
void free_clusters(struct cluster_t *clusters, int ncl, int freeclarr)
{
    for (int i = 0; i < ncl; i++) {
        free(clusters[i].obj);
    }
    if (freeclarr)
        free(clusters);
}

/*
 Ze souboru 'filename' nacte objekty. Pro kazdy objekt vytvori shluk a ulozi
 jej do pole shluku. Alokuje prostor pro pole vsech shluku a ukazatel na prvni
 polozku pole (ukalazatel na prvni shluk v alokovanem poli) ulozi do pameti,
 kam se odkazuje parametr 'arr'. Funkce vraci pocet nactenych objektu (shluku).
 V pripade nejake chyby uklada do pameti, kam se odkazuje 'arr', hodnotu NULL.
*/
int load_clusters(char *filename, struct cluster_t **arr)
{
    assert(arr != NULL);

    int count = 0;
    
    FILE *file = fopen(filename, "r");

    checkforerror(file == NULL, "error: file not found\n");
    checkforerror(fscanf(file, "count=%d\n", &count) != 1, "error: invalid count\n");
    checkforerror(count < 1, "error: count can't be below one\n");

    struct cluster_t poleshluku[count];
    int ids[count];

    int id = 0; float x = 0, y = 0, temp_id = 0;
    for (int i = 0; i < count; i++) {
        struct cluster_t temp_cluster;
        if (fscanf(file, "%f %f %f\n", &temp_id, &x, &y) != 3) { // if there is no 3 arguments
            fprintf(stderr, "error: invalid input (less than 3 numbers in line %d)\n", i+2);
            free_clusters(poleshluku, i, 0);
            fclose(file);
            exit(1);
        }
        if (trunc(temp_id) != temp_id || trunc(x) != x || trunc(y) != y) { // if id is not integer
            fprintf(stderr, "error: data are not integer (on line %d)\n", i+2);
            free_clusters(poleshluku, i, 0);
            fclose(file);
            exit(1);
        }
        id = temp_id;
        if (id < 0 || x > 1000 || y > 1000 || x < 0 || y < 0) {
            fprintf(stderr, "error: invalid input (on line %d)\n", i+2);
            free_clusters(poleshluku, i, 0);
            fclose(file);
            exit(1);
        }
        for (int j = 0; j < i; j++) { // test for unique id
            if (ids[j] == id) {
                fprintf(stderr, "error: id already exists\n");
                free_clusters(poleshluku, i, 0);
                fclose(file);
                exit(1);
            }
        }
        struct obj_t temp_object = {.id=id, .x=x, .y=y};
        init_cluster(&temp_cluster, 1);
        append_cluster(&temp_cluster, temp_object);
        poleshluku[i] = temp_cluster;
        ids[i] = id;
    }

    struct cluster_t *clusterarray = malloc(sizeof(poleshluku));
    for (int i = 0; i < count; i++) {
        clusterarray[i] = poleshluku[i];
    }
    *arr = clusterarray;

    fclose(file);
    return count;
}

/*
 Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
 Tiskne se prvnich 'narr' shluku.
*/
void print_clusters(struct cluster_t *carr, int narr)
{
    printf("Clusters:\n");
    for (int i = 0; i < narr; i++)
    {
        printf("cluster %d: ", i);
        print_cluster(&carr[i]);
    }
}

/*
 Pocita vzdalenost dvou shluku,
 pomoci metody nejvzdálenějšího souseda.
*/
float cd_maximum_objd(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float dist = obj_distance(&c1->obj[0], &c2->obj[0]);

    for (int i = 0; i < c1->size; i++) {
        for (int j = 0; j < c2->size; j++) {
            dist = (dist < obj_distance(&c1->obj[i], &c2->obj[j])) ? obj_distance(&c1->obj[i], &c2->obj[j]) : dist;
        }
    }
    return dist;
}

// k-means clustering

/*
  Extracts all the objects from the cluster array
  and puts them into one array.
*/
void extract_objects(struct cluster_t *arr, int narr, struct obj_t *objarr) {
    assert(arr != NULL);
    assert(objarr != NULL);
    assert(narr > 0);
    for (int i = 0; i < narr; i++) {
        objarr[i] = arr[i].obj[0];
    }
}

/*
Compares two obj_t 
*/
int compare_objects(struct obj_t *a, struct obj_t *b, int k) {
    assert(a != NULL);
    assert(b != NULL);
    assert(k > 0);
    for (int i = 0; i < k; i++) {
        if (a[i].x == b[i].x && a[i].y == b[i].y) return 1;
    }
    return 0;
}

float squared_obj_distance(struct obj_t *o1, struct obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    return pow((o1->x - o2->x), 2)+pow((o1->y - o2->y), 2);
}

int nearest_point(struct obj_t *givenobj, struct obj_t *centroids, int k) {
    assert(k > 1);
    float dist = squared_obj_distance(givenobj, &centroids[0]);
    int idnearestcl = 0;
    for (int j = 1; j < k; j++)  {
        if (squared_obj_distance(givenobj, &centroids[j]) < dist) {
            dist = squared_obj_distance(givenobj, &centroids[j]);
            idnearestcl = j;
        }
    }
    return idnearestcl;
}

struct obj_t *k_meanspp(struct obj_t *objarr, int narr, int k) { // k-means++
    assert(objarr != NULL);
    assert(narr > 0);
    assert(k > 0);

    struct obj_t *centroids = malloc(sizeof(struct obj_t)*k);
    centroids[0] = objarr[rand() % narr];
    for (int i = 1; i < k; i++) {
        float sumofdist = 0;
        float *distances = malloc(sizeof(float)*narr);
        for (int j = 0; j < narr; j++) {
            float mindist = squared_obj_distance(&objarr[j], &centroids[0]);
            for (int l = 1; l < i; l++) {
                if (squared_obj_distance(&objarr[j], &centroids[l]) < mindist) {
                    mindist = squared_obj_distance(&objarr[j], &centroids[l]);
                }
            }
            distances[j] = mindist;
            sumofdist += mindist;
        }
        float random = ((float)rand()/RAND_MAX)*sumofdist;
        float sum = 0;
        for (int j = 0; j < narr; j++) {
            sum += distances[j];
            if (sum >= random) {
                centroids[i] = objarr[j];
                break;
            }
        }
        free(distances);
    }
    return centroids;
}

struct obj_t *k_meanspp_roll(struct obj_t *objarr, int narr, int k, int acc) {
    struct obj_t *centroids = k_meanspp(objarr, narr, k);
    float mean = 0.0;
    for (int j = 0; j < narr; j++) {
            mean += squared_obj_distance(&objarr[j], &centroids[nearest_point(&objarr[j], centroids, k)]);
    }
    for (int i = 1; i < acc; i++) {
        struct obj_t *tempcentroids = k_meanspp(objarr, narr, k);
        float tempmean = 0.0;
        for (int j = 0; j < narr; j++) {
            tempmean += squared_obj_distance(&objarr[j], &tempcentroids[nearest_point(&objarr[j], tempcentroids, k)]);
        }
        if (tempmean < mean) {
            mean = tempmean;
            centroids = tempcentroids;
        }
    }
    // fprintf(stderr, "inert: %f\n", mean);
    return centroids;
}

struct cluster_t *k_means(struct obj_t *objarr, int narr, int k, int acc) {
    struct obj_t *centroids = k_meanspp_roll(objarr, narr, k, acc);
    int iterations = 0; int numofbestcl;
    struct obj_t *oldcentroids = calloc(k, sizeof(struct obj_t));
    struct cluster_t *labels = calloc(k, sizeof(struct cluster_t));
    for (int i = 0; i < k; i++) init_cluster(&labels[i], 1);
    while (!compare_objects(centroids, oldcentroids, k) && iterations < 200) {
        for (int i = 0; i < k; i++) clear_cluster(&labels[i]);
        for (int i = 0; i < narr; i++) {
            numofbestcl = nearest_point(&objarr[i], centroids, k);
            append_cluster(&labels[numofbestcl], objarr[i]);
        }

        memcpy(oldcentroids, centroids, sizeof(struct obj_t) * k);
        struct obj_t newcentroid;
        for (int i = 0; i < k; i++) {
            newcentroid = labels[i].obj[0];
            for (int j = 1; j < labels[i].size; j++) {
                newcentroid.x += labels[i].obj[j].x;
                newcentroid.y += labels[i].obj[j].y;
            }
            newcentroid.x /= labels[i].size;
            newcentroid.y /= labels[i].size;
            centroids[i] = newcentroid;
        }
        iterations++;
    }
    free(centroids);
    free(oldcentroids);
    return labels;
}

int main(int argc, char *argv[])
{
    struct cluster_t *clusters;

    int clusterstosortnum;
    float (*dist_fun)(struct cluster_t *, struct cluster_t *);
    if (argc == 2) {
        clusterstosortnum = 1;
        dist_fun = cluster_distance;
    }
    else if (argc == 3) {
        checkforerror(atoi(argv[2]) == 0, "invalid argument (not number)\n")
        checkforerror(atoi(argv[2]) != atof(argv[2]), "invalid argument (must be integer)\n")
        dist_fun = cluster_distance;
        clusterstosortnum = atoi(argv[2]);
    }
    else if (argc == 4 && strcmp(argv[1], "-k") == 0) { // k-means
        clusterstosortnum = atoi(argv[3]);
        int count = load_clusters(argv[2], &clusters);
        srand(time(NULL));
        struct obj_t *objarr = malloc(sizeof(struct obj_t) * count);
        extract_objects(clusters, count, objarr);
        struct cluster_t *sortedarr = k_means(objarr, count, clusterstosortnum, 1000);
        print_clusters(sortedarr, clusterstosortnum);
        for (int i = 0; i < clusterstosortnum; i++) {
            free(clusters[i].obj);
        }
        free(clusters);
        return 0;
    }
    else if (argc == 4 && strcmp(argv[3], "-nvzd") == 0) { // nejvzdalen. soused
        dist_fun = &cd_maximum_objd;
        clusterstosortnum = atoi(argv[2]);
    }
    else {
        fprintf(stderr, "incorrect number of arguments! try: ./cluster 'file' [clustersnum]\n");
        return 1;
    }

    int c1, c2;

    int count = load_clusters(argv[1], &clusters);
    if (count == -1) 
        return -1;
    int sizeofarr = count;
    for (int i = 0; i < count - clusterstosortnum; i++) {
        find_neighbours(clusters, sizeofarr, &c1, &c2, dist_fun);
        merge_clusters(&clusters[c1], &clusters[c2]);
        sizeofarr = remove_cluster(clusters, sizeofarr, c2);
    }
    for (int i = 0; i < clusterstosortnum; i++) {
        sort_cluster(&clusters[i]);
    }
    print_clusters(clusters, sizeofarr);

    free_clusters(clusters, sizeofarr, 1);
}