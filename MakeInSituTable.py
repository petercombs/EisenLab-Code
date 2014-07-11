from __future__ import print_function
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl
from tkinter import Tk,  Frame, Button, IntVar, Label, Menu
import tkinter.filedialog
from PIL import Image, ImageTk
import urllib.request, urllib.parse, urllib.error, time
from os import path, makedirs
from math import ceil, floor, sqrt
from idlelib.ToolTip import ToolTip
from sys import stdout

if not path.isfile('figures/insitu_images.csv'):
    print("Getting necessary files...")
    from subprocess import call
    try:
        makedirs('figures/BDGP')
    except OSError:
        pass
    urllib.request.urlretrieve('http://insitu.fruitfly.org/insitu-mysql-dump/insitu_images.csv.gz',
                       'figures/insitu_images.csv.gz')
    call(['gunzip', 'figures/insitu_images.csv.gz'])


BDGP_names = pd.read_csv('figures/insitu_images.csv', header=None)
BDGP_names.columns =  ['Gene name', 'Name2', 'Name3', 'FBgn',
                                  'SomethingElse', 'imgfile', 'stage',
                                   'view', 'ap_invert', 'hilo']

st_5_with_lat = BDGP_names[(BDGP_names['stage'] == 2)]
has_img = set(st_5_with_lat['Name2']).union(set(st_5_with_lat['FBgn']))

def imget(imname):
    """ Use cached, or fetch an image from FlyExpress

    Assumes that the image name is one from BDGP, in which case
    it's pretty easy to look at the source of the FlyExpress
    report pages and see what the format is.
    """
    im_basename = path.splitext(path.basename(imname))[0]
    filename = path.join('figures', 'BDGP', im_basename+'.bmp')
    if not path.exists(filename):
        base_web = "http://www.flyexpress.net/ZOOX4_DBImages/BDGP/thumbnails/%s_s.bmp"
        print("1 second delay to avoid spamming server")
        time.sleep(1)
        ret = urllib.request.urlretrieve(base_web % im_basename, filename)
    return (filename)

def imget_bdgp(imname):
    """ Use cached, or fetch an image from FlyExpress

    Assumes that the image name is one from BDGP, in which case
    it's pretty easy to look at the source of the FlyExpress
    report pages and see what the format is.
    """
    im_basename = path.splitext(path.basename(imname))[0]
    filename = path.join('figures', 'BDGP', im_basename+'.jpg')
    if not path.exists(filename):
        base_web = "http://insitu.fruitfly.org/insitu_image_storage/thumbnails/%s"
        print("1 second delay to avoid spamming server")
        time.sleep(1)
        ret = urllib.request.urlretrieve(base_web % imname, filename)

    return (filename)
def adjust_rows_cols(n, rows, cols):
    if rows > 0 and cols > 0:
        return rows, cols
    elif rows > 0:
        return rows, int(ceil(n / float(rows)))
    elif cols > 0:
        return int(ceil(n / float(cols))), cols
    else:
        return ceil(sqrt(n)), floor(sqrt(n))

class App:
    def __init__(self, master, cluster_fname=None, fig_cols = 4, fig_rows = None):
        if not cluster_fname:
            cluster_fname = tkinter.filedialog.askopenfilename(parent=master,
                                                         title="List of genes",
                                                         initialdir=path.abspath('./'))
            print(cluster_fname)
        master.wm_title(cluster_fname)
        self.all_in_cluster = [line.strip() for line in open(cluster_fname)]
        print("No images", [item for item in self.all_in_cluster
                            if item not in has_img])
        self.imgs_in_cluster = [item for item in self.all_in_cluster
                                if item in has_img]
        print(self.imgs_in_cluster)
        print(len(self.imgs_in_cluster),  len(self.all_in_cluster), end=' ')
        print(len(self.imgs_in_cluster) / float(len(self.all_in_cluster)))
        self.img_lists = [[x for x in
                           st_5_with_lat[(st_5_with_lat['Name2'] == gene)
                                        +(st_5_with_lat['FBgn'] == gene)]['imgfile']]
                          for gene in self.imgs_in_cluster]
        num_genes = len(self.imgs_in_cluster)
        self.fig_rows, self.fig_cols = adjust_rows_cols(num_genes,
                                                        fig_rows or 0, fig_cols or 0)
        self.outfn = path.basename(cluster_fname)[:-4]+'.png'


        buttons = Frame(master)
        buttons.grid(row=1, column=1)
        savebutton = Button(buttons, text="Save", command=self.save)
        savebutton.grid(row=2, column = 2)
        Button(buttons, text="Add Row", command=self.add_row).grid(row=1,column=3)
        Button(buttons, text="Add Col", command=self.add_col).grid(row=3,column=3)
        Button(buttons, text="Remove Row",
               command=self.remove_row).grid(row=1,column=1)
        Button(buttons, text="Remove Col",
               command=self.remove_col).grid(row=3,column=1)

        images = Frame(master)
        images.grid(row=2, column=1)


        self.button_selected = IntVar(master, -1)
        self.images_selected = [0 for list in self.img_lists]
        out_name = '{}.png.txt'.format(path.splitext(path.basename(cluster_fname))[0])
        print(out_name)
        if path.exists(out_name):
            for i, line in enumerate(open(out_name)):
                self.images_selected[i] = int(line.split()[0].strip())
        self.radio_buttons = []
        self.image_files = []
        self.x_flipped = []
        for i, gene in enumerate(self.imgs_in_cluster):
            for image_file_name in self.img_lists[i]:
                try:
                    raise IOError
                    image_file = imget(image_file_name)
                    import_image = Image.open(image_file)
                except (IOError, OSError):
                    try:
                        image_file = imget_bdgp(image_file_name)
                        import_image = Image.open(image_file)
                    except Exception as exc:
                        print(exc)
                        print(image_file_name, gene)
                        raise exc
                    import_image = import_image.resize((180,80), Image.ANTIALIAS)
                except urllib.error.HTTPError as exc:
                    print(image_file_name, gene)
                    raise exc
                embryo_image = ImageTk.PhotoImage(import_image)

            label = Label(image=embryo_image)
            label.image = embryo_image
            b = Button(images,
                       image = label.image,
                       width = 200,
                       height= 100,
                       text=gene,
                       #text = image_file,
                       #variable=self.button_selected,
                       #command=self.handle_button_click,
                       #value = i+1,
                       #indicatoron=0,
                      )
            #ToolTip(b, gene)
            b.bind("<Button-1>", self.handle_button_click)
            b.bind("<Button-2>", self.handle_button_click)
            b.grid(row = i // self.fig_cols, column = i % self.fig_cols)
            #w = Tix.Balloon(master, message=gene)
            #w.bind_widget(b)

            self.image_files.append(image_file)
            self.radio_buttons.append(b)
            self.x_flipped.append(False)


    def save(self):
        outfh = open(self.outfn+'.txt', 'w')
        [outfh.write(str(val) + '\n') for val in self.images_selected]
        outfh.close()
        fig = mpl.figure(figsize=(1.4*self.fig_cols, self.fig_rows))
        for i, gene in enumerate(self.image_files):
            ax = fig.add_subplot(self.fig_rows, self.fig_cols, i+1)
            try:
                ax.imshow(mpl.imread(gene), origin='lower')
                print(gene)
                ax.set_title(self.imgs_in_cluster[i])
            except IOError:
                print(i // int(self.fig_cols), i % int(self.fig_cols), "Can't find image for", gene)
            mpl.xticks([], [])
            mpl.yticks([], [])
            if self.x_flipped[i]:
                ax.invert_xaxis()
            for spine in ax.spines.values():
                spine.set_visible(False)
                pass
        fig.tight_layout()
        mpl.savefig(path.join('figures', self.outfn), dpi=300)
        mpl.savefig(path.join('figures', self.outfn + '.eps'), dpi=300)

    def handle_button_click(self, event):

        val = self.button_selected.get() - 1
        val = self.radio_buttons.index(event.widget)
        if event.state & 0x001:  # Shift key down
            self.x_flipped[val] = not self.x_flipped[val]
            try:
                raise IOError
                image_file = imget(self.img_lists[val][self.images_selected[val]])
                import_image = Image.open(image_file)
            except IOError:
                image_file = imget_bdgp(self.img_lists[val][self.images_selected[val]])
                import_image = Image.open(image_file)
                import_image = import_image.resize((180,80), Image.ANTIALIAS)
            import_image = import_image.rotate(180*self.x_flipped[val])
            embryo_image = ImageTk.PhotoImage(import_image)
            label = Label(image=embryo_image)
            label.image = embryo_image
            self.radio_buttons[val].configure(
                #text=image_file,
                image=label.image,
            )
        elif event.state & 0x004: # Control key down
            popup = Menu(root, tearoff=0)
            popup.add_command(label=self.radio_buttons[val].cget('text'))
            popup.add_separator()
            popup.add_command(label="Invert X")
            popup.add_command(label="Invert Y")
            try:
                popup.tk_popup(event.x_root, event.y_root, 0)
            finally:
                # make sure to release the grab (Tk 8.0a1 only)
                popup.grab_release()

        elif val >= 0:
            self.images_selected[val] = ((self.images_selected[val] + 2*event.num
                                         - 3) %
                                         len(self.img_lists[val]))
            self.x_flipped[val] = False
            if self.images_selected[val] == 0:
                stdout.write('\a')
                stdout.flush()

            try:
                raise IOError
                image_file = imget(self.img_lists[val][self.images_selected[val]])
                import_image = Image.open(image_file)
            except IOError:
                image_file = imget_bdgp(self.img_lists[val][self.images_selected[val]])
                import_image = Image.open(image_file)
                import_image = import_image.resize((180,80), Image.ANTIALIAS)
            embryo_image = ImageTk.PhotoImage(import_image)
            label = Label(image=embryo_image)
            label.image = embryo_image
            self.radio_buttons[val].configure(
                #text=image_file,
                image=label.image,
            )
            self.image_files[val] = image_file
        self.button_selected.set(-1)

    def add_row(self):
        self.fig_rows += 1
        self.fig_rows, self.fig_cols = adjust_rows_cols(len(self.imgs_in_cluster),
                                                        self.fig_rows, 0)
        print(self.fig_rows, self.fig_cols)
        self.redraw_buttons()

    def remove_row(self):
        self.fig_rows = max(self.fig_rows - 1, 1)
        self.fig_rows, self.fig_cols = adjust_rows_cols(len(self.imgs_in_cluster),
                                                        self.fig_rows, 0)
        print(self.fig_rows, self.fig_cols)
        self.redraw_buttons()

    def add_col(self):
        self.fig_cols += 1
        self.fig_rows, self.fig_cols = adjust_rows_cols(len(self.imgs_in_cluster),
                                                        0, self.fig_cols)
        print(self.fig_rows, self.fig_cols)
        self.redraw_buttons()

    def remove_col(self):
        self.fig_cols = max(self.fig_cols - 1, 1)
        self.fig_rows, self.fig_cols = adjust_rows_cols(len(self.imgs_in_cluster),
                                                        0, self.fig_cols)
        print(self.fig_rows, self.fig_cols)
        self.redraw_buttons()

    def redraw_buttons(self):
        for i, b in enumerate(self.radio_buttons):
            b.grid(row = i // self.fig_cols, column = i % self.fig_cols)



root = Tk()
#root.tk.call("load", "", "Tix")
#app = App(root, 'Cluster12.txt')
app = App(root)
root.mainloop()


