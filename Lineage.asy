import roundedpath;
import patterns;
import Distributions;

struct mutation{
    static void draw(path p, int marker, pen color = red){
        draw(midpoint(p), scale(10) * unitcircle, linewidth(5bp) + color);
        fill(midpoint(p), scale(10) * unitcircle, gray(1) + color);
        label(scale(3) * string(marker + 1),
              shift((0.8, -0.8)) * midpoint(p),
              p = AvantGarde());
    }
}

struct Marker {
    bool3 status; // true -> derived, false -> wild, default -> non-ancestral.
    pair position;

    // Box representing the marker.
    pen box_pen;

    // Colors.
    pen derived_color;
    pen wild_color;
    pen nonancestral_color;

    bool3 getstatus() {
        return this.status;
    }

    void setnonancestral() {
        this.status = default;
    }

    void shiftby(pair p) {
        this.position += p;
    }

    void operator init(bool3 status,
                       pair pos,
                       pen derived_color = red,
                       pen wild_color = blue,
                       pen nonancestral_color = lightblue){
        this.status = status;
        this.position = pos;

        this.box_pen = linewidth(3bp);
        if (status == default)
            box_pen += nonancestral_color;
        else
            box_pen += status ? derived_color : wild_color;

        this.derived_color = derived_color;
        this.wild_color = wild_color;
        this.nonancestral_color = nonancestral_color;
    }

    void draw(){
        path box = shift(this.position) * roundedpath(polygon(4), 0.2);

        // Hide the line in the background.
        fill(box, white);

        // Pattern for non ancestral material.
        add("hatch", hatch(p = linewidth(4bp) + nonancestral_color));

        if (this.status == default) {
            fill(box, pattern("hatch"));
        }
        else {
            fill(box, gray(1) + (this.status ? derived_color : wild_color));

            string marker_text = this.status ? "1" : "0";

            pen marker_pen = AvantGarde(series = "b") + fontsize(30pt) +
                (this.status ? derived_color : wild_color);

            label(marker_text, this.position, p = marker_pen);
        }

        draw(box, this.box_pen);
    }
}

struct Sequence {
    Marker[] markers;
    Sequence[] _children;

    // Recombination stuff.
    int _recombination_point;

    // Mutation stuff.
    string _mutation_branch;
    int _mutation_locus;

    pair _rightpoint;
    pair _north_anchor;
    pair _south_anchor;
    pair _position;

    // Colors.
    pen derived_color;
    pen wild_color;
    pen nonancestral_color;

    int nbmarkers() {
        return this.markers.length;
    }

    bool3 getstatus(int marker) {
        return this.markers[marker].getstatus();
    }

    void showguide(string s) {
        label(scale(7.5) * s,
              midpoint(this._north_anchor -- this._south_anchor));
    }

    bool3[] status() {
        bool3[] ret = array(this.markers.length, default);
        for (int k = 0; k < ret.length; ++k) {
            bool3 status = this.markers[k].getstatus();
            if (status != default)
                ret[k] = status;
        }

        return ret;
    }

    // Base constructor.
    void operator init(bool3[] status,
                       pair pos,
                       pen derived_color = red,
                       pen wild_color = blue,
                       pen nonancestral_color = lightblue,
                       Sequence[] children = {},
                       string mutation_branch = "",
                       int mutation_locus = -1,
                       int recombination_point = -1){
        this.derived_color = derived_color;
        this.wild_color = wild_color;
        this.nonancestral_color = nonancestral_color;

        this._children = children;

        // Recombination stuff.
        this._recombination_point = recombination_point;

        // Mutation stuff.
        this._mutation_branch = mutation_branch;
        this._mutation_locus = mutation_locus;

        this._position = pos;

        this._rightpoint =
            shift(_position) * ((status.length - 1) * (sqrt(2) + 0.25), 0);

        this._north_anchor =
            shift(_position) *
            (this._rightpoint.x, this._rightpoint.y + (sqrt(2) + 0.5)) / 2;
        this._south_anchor =
            shift(_position) *
            (this._rightpoint.x, this._rightpoint.y - (sqrt(2) + 0.5)) / 2;
        
        for (int kk = 0; kk < status.length; ++kk){
            pair final_position = shift(_position) * (kk * (sqrt(2) + 0.25), 0);
            this.markers.push(Marker(status[kk],
                                     final_position,
                                     derived_color,
                                     wild_color,
                                     nonancestral_color));
        }
    }

    // Deep copy constructor.
    void operator init(Sequence other){
        this.markers = copy(other.markers);
        this._children = copy(other._children);
        this._recombination_point = other._recombination_point;
        this._mutation_branch = other._mutation_branch;
        this._mutation_locus = other._mutation_locus;
        this._rightpoint = other._rightpoint;
        this._north_anchor = other._north_anchor;
        this._south_anchor = other._south_anchor;
        this._position = other._position;
    }

    // Get the width of a sequence.
    real width(){
        return this._rightpoint.x - this._position.x;
    }

    static Sequence right_of(Sequence other,
                             bool3[] status,
                             pen derived_color = red,
                             pen wild_color = blue,
                             pen nonancestral_color = lightblue,
                             real distance = 5){

        pair pos = shift(distance * sqrt(2) / 2, 0) * other._rightpoint;
        return Sequence(status,
                        pos,
                        derived_color,
                        wild_color,
                        nonancestral_color);
    }

    void draw(bool mutations = true,
              pen edges_color = grey,
              pen derived_color = red){
        // Draw sequences.
        draw(shift(_position) * (0, 0) -- _rightpoint,
             linewidth(3bp) + edges_color);

        for (Marker m : this.markers)
            m.draw();

        // Draw edges.
        path[] edges;
        for (Sequence child : _children)
            edges.push(child._north_anchor -- this._south_anchor);

        for (path edge : edges)
            draw(edge, linewidth(4bp) + edges_color);

        // Draw mutations.
        if (mutations){
            if (this._mutation_branch == "left")
                mutation.draw(edges[0], this._mutation_locus, derived_color);
            else if (this._mutation_branch == "right")
                mutation.draw(edges[1], this._mutation_locus, derived_color);
        }
    }

    static void draw(Sequence[] sequences,
                     bool mutations = true,
                     pen edges_color = grey,
                     pen derived_color = red){
        for (Sequence sequence : sequences)
            sequence.draw(mutations, edges_color, derived_color);
    }

    Sequence coalesce(Sequence other, real length = 3 * sqrt(2)){
        // Determine if this is at the left of other.
        bool this_left = this._position.x < other._position.x;

        // Compute the new sequence.
        bool3[] new_status = new bool3[this.nbmarkers()];
        for (int kk = 0; kk < this.nbmarkers(); ++kk){
            bool3 this_status = this.getstatus(kk);
            bool3 other_status = other.getstatus(kk);

            if (this_status != default && other_status != default)
                new_status[kk] = this_status && other_status;
            else if (this_status != default)
                new_status[kk] = this_status;
            else if (other_status != default)
                new_status[kk] = other_status;
            else
                new_status[kk] = default;
        }

        // Determine if a mutation took place.
        int mutation_locus = -1;
        for (int kk = 0; kk < this.nbmarkers(); ++kk){
            // If one of the coalescing vertices is not ancestral for marker of
            // interest, no mutation could have occured.
            if (this.getstatus(kk) == default || other.getstatus(kk) == default)
                continue;

            // TODO: do not assume infinite sites model.
            if (this.getstatus(kk) ^ other.getstatus(kk)){
                mutation_locus = kk;
                break;
            }
        }

        string mutation_branch = "";
        if (mutation_locus >= 0) {
            if (this.getstatus(mutation_locus))
                mutation_branch = "left";
            else
                mutation_branch = "right";
        }

        // Compute the position of the new sequence.
        real new_x = (this._position.x + other._position.x) / 2;
        real new_y = max(this._position.y, other._position.y);
        pair new_position =
            shift(0, length) * (new_x, new_y);

        // Children of the new node.
        Sequence[] new_children = {this, other};

        return Sequence(new_status,
                        new_position,
                        this.derived_color,
                        this.wild_color,
                        this.nonancestral_color,
                        children = new_children,
                        mutation_branch = mutation_branch,
                        mutation_locus = mutation_locus);
    }

    Sequence mutate(int mutating_marker, int n, real length = 3 * sqrt(2)){
        // Compute the new sequence.
        bool3[] new_status = new bool3[this.nbmarkers()];
        for (int k = 0; k < this.nbmarkers(); ++k)
            new_status[k] = this.getstatus(k);

        new_status[mutating_marker] = false;

        // Compute the position of the new sequence.
        pair new_position =
            shift(0, length) * this._position;

        // Child of the new node.
        Sequence[] new_child = {this};

        return Sequence(new_status,
                        new_position,
                        children = new_child,
                        mutation_branch = "left",
                        mutation_locus = mutating_marker);
    }

    // Recombination point is left inclusive!!!
    Sequence[] recombine(int breakpoint,
                         real length = 3 * sqrt(2)){
        // Compute new status.
        bool3[][] newstatus = {this.status(), this.status()};
        newstatus[0][breakpoint + 1:] =
            array(newstatus[0].length - breakpoint - 1, default);
        newstatus[1][:breakpoint + 1] =
            array(breakpoint + 1, default);

        // Compute new positions.
        real xshift = this.width();
        pair[] newpositions = {shift(-xshift, length) * this._position,
                               shift(xshift, length) * this._position};


        Sequence[] ret = new Sequence[2];
        for (int k = 0; k < 2; ++k)
            ret[k] = Sequence(newstatus[k],
                              newpositions[k],
                              this.derived_color,
                              this.wild_color,
                              this.nonancestral_color,
                              new Sequence[]{this});

        return ret;
    }
}

struct Arg {
    Sequence[] sequences;

    // Should mutation edges be marked?
    bool draw_mutations;

    // Colors.
    pen edges_color;
    pen derived_color;
    pen wild_color;
    pen nonancestral_color;

 void operator init(bool3[] leaf,
                    bool draw_mutations = true,
                    pen edges_color = grey,
                    pen derived_color = red,
                    pen wild_color = blue,
                    pen nonancestral_color = lightblue) {
        this.draw_mutations = draw_mutations;
        this.edges_color = edges_color;
        this.derived_color = derived_color;
        this.wild_color = wild_color;
        this.nonancestral_color = nonancestral_color;

        this.sequences.push(Sequence(leaf,
                                     (0, 0),
                                     derived_color,
                                     wild_color,
                                     nonancestral_color));
    }

    void newleaf(bool3[] status, int n = 1) {
        for (int k = 0; k < n; ++k) {
            Sequence last_sequence = this.sequences[this.sequences.length - 1];
            Sequence sequence = Sequence.right_of(last_sequence,
                                                  status,
                                                  this.derived_color,
                                                  this.wild_color,
                                                  this.nonancestral_color);
            this.sequences.push(sequence);
        }
    }

    void coalesce(int seq_idx1, int seq_idx2) {
        Sequence coal_seq = this.sequences[seq_idx1].coalesce(this.sequences[seq_idx2]);
        this.sequences.push(coal_seq);
    }

    void recombine(int seq_idx, int breakpoint) {
        this.sequences.append(this.sequences[seq_idx].recombine(breakpoint));
    }

    void draw(bool guides = false) {
        Sequence.draw(sequences,
                      this.draw_mutations,
                      this.edges_color,
                      this.derived_color);

        if (guides)
            for (int k = 0; k < this.sequences.length; ++k) {
                sequences[k].showguide(string(k));
            }
    }
}
