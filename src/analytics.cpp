#include <cstdlib>
#include <cstring>
#include <clocale>
#include "spectrum.h"
#include "plots.h"
#include <gtk/gtk.h>

GtkWidget *A, *P, *T, *K, *Ecut, *Elem;

static void
start (GtkWidget *widget,
       gpointer   data)
{
    bool a, p, t, k;
    a = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(A));
    p = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(P));
    t = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(T));
    k = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(K));

    int N = 1000;
    int M = 5000;


    double l;
    int z;
    const gchar *el = gtk_entry_get_text(GTK_ENTRY(Elem));
    if (!strcmp(el, "Si"))
    {
        z = 14;
        l = 1e-4;
    }
    else if (!strcmp(el, "Ge"))
    {
        z = 32;
        l = 6e-5;
    }
    else
    {
        gtk_entry_set_text(GTK_ENTRY(Elem), "Si or Ge");
        gtk_widget_grab_focus(Elem);
        return;
    }

    double ecut = atof(gtk_entry_get_text(GTK_ENTRY(Ecut)));
    if (a)
        analytical(z, M, N, l, ecut);

    if (t && z == 32)
        table(z, M, N, l, ecut);

    if (p)
        approximation(z, M, N, l, ecut);

    plot_analytics(a, p, t && z == 32);

    if (k)
    {
        quit_function(z, M, N, l, ecut);
        plot_k();
    }
}

int
main (int   argc,
      char *argv[])
{
    GtkWidget *window;
    GtkWidget *grid;
    GtkWidget *button;
    GtkWidget *label;

    /* This is called in all GTK applications. Arguments are parsed
    * from the command line and are returned to the application.
    */
    gtk_init (&argc, &argv);
    setlocale(LC_ALL, "en_US.UTF-8");

    /* create a new window, and set its title */
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title (GTK_WINDOW (window), "Analytics");
    g_signal_connect (window, "destroy", G_CALLBACK (gtk_main_quit), NULL);
    gtk_container_set_border_width (GTK_CONTAINER (window), 20);

    /* Here we construct the container that is going pack our buttons */
    grid = gtk_grid_new ();

    /* Pack the container in the window */
    gtk_container_add (GTK_CONTAINER (window), grid);

    label = gtk_label_new("Emitter:");
    gtk_grid_attach (GTK_GRID (grid), label, 0, 0, 1, 1);

    Elem = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(Elem), "Si or Ge");
    gtk_grid_attach (GTK_GRID (grid), Elem, 0, 1, 1, 1);

    label = gtk_label_new("Ecut, eV:");
    gtk_grid_attach (GTK_GRID (grid), label, 0, 2, 1, 1);

    Ecut = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(Ecut), "500");
    gtk_grid_attach (GTK_GRID (grid), Ecut, 0, 3, 1, 1);

    label = gtk_label_new("Types of data:");
    gtk_grid_attach (GTK_GRID (grid), label, 0, 4, 1, 1);

    A = gtk_check_button_new_with_label("analytical (A)");
    gtk_grid_attach (GTK_GRID (grid), A, 0, 5, 1, 1);

    T = gtk_check_button_new_with_label("table (T)");
    gtk_grid_attach (GTK_GRID (grid), T, 0, 6, 1, 1);

    P = gtk_check_button_new_with_label("approximation (P)");
    gtk_grid_attach (GTK_GRID (grid), P, 0, 7, 1, 1);

    label = gtk_label_new("Another options:");
    gtk_grid_attach (GTK_GRID (grid), label, 0, 8, 1, 1);
    K = gtk_check_button_new_with_label("Exit function K(t)");
    gtk_grid_attach (GTK_GRID (grid), K, 0, 9, 1, 1);

    button = gtk_button_new_with_label ("Start!");
    g_signal_connect (button, "clicked", G_CALLBACK (start), NULL);
    gtk_grid_attach (GTK_GRID (grid), button, 0, 10, 1, 1);

    /* Now that we are done packing our widgets, we show them all
    * in one go, by calling gtk_widget_show_all() on the window.
    * This call recursively calls gtk_widget_show() on all widgets
    * that are contained in the window, directly or indirectly.
    */
    gtk_widget_show_all (window);

    /* All GTK applications must have a gtk_main(). Control ends here
    * and waits for an event to occur (like a key press or a mouse event),
    * until gtk_main_quit() is called.
    */
    gtk_main ();

    return 0;
}
