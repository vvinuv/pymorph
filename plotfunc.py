import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


class FitsPlotter:

    def __init__(self, output_file, mask_file):
        self.output_file = output_file
        self.mask_file = mask_file

        self.image, self.model, self.residual = self.load_fits()
        self.mask = fits.open(mask_file)[0].data

        self.image_m, self.valid = self.apply_mask(self.image)
        self.model_m, _ = self.apply_mask(self.model)
        self.residual_m, _ = self.apply_mask(self.residual)

    # -------- Load FITS --------
    def load_fits(self):
        hdul = fits.open(self.output_file)
        return hdul[1].data, hdul[2].data, hdul[3].data

    # -------- Apply Mask --------
    def apply_mask(self, data):
        valid = (self.mask == 0)
        return np.where(valid, data, np.nan), valid

    # -------- Plot Image --------
    def plot_image(self, ax, data, title, cmap='gray'):
        ny, nx = data.shape
        x = np.arange(nx) - nx // 2
        y = np.arange(ny) - ny // 2

        vmin = np.nanpercentile(data, 1)
        vmax = np.nanpercentile(data, 99)

        ax.imshow(data, origin='lower',
                  extent=[x.min(), x.max(), y.min(), y.max()],
                  cmap=cmap, vmin=vmin, vmax=vmax)

        ax.set_title(title, size=15)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_aspect('equal')

    # -------- Residual --------
    def plot_residual(self, ax, residual, title):
        ny, nx = residual.shape
        x = np.arange(nx) - nx // 2
        y = np.arange(ny) - ny // 2

        vmax = np.nanpercentile(np.abs(residual), 99)

        ax.imshow(residual, origin='lower',
                  extent=[x.min(), x.max(), y.min(), y.max()],
                  cmap='coolwarm', vmin=-vmax, vmax=vmax)

        ax.set_title(title, size=15)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)

        ax.set_aspect('equal')

    # -------- Radial Profile --------
    def radial_profile(self, center=None):
        image = self.image_m
        model = self.model_m
        valid = self.valid

        y, x = np.indices(image.shape)

        if center is None:
            center = np.array(image.shape) // 2

        r = np.sqrt((x - center[1])**2 + (y - center[0])**2).astype(int)

        r_flat = r[valid]
        img_flat = image[valid]
        mod_flat = model[valid]

        r_max = r_flat.max()
        bins = np.arange(1, r_max)

        img_prof, mod_prof, img_err = [], [], []

        for i in bins:
            mask_bin = (r_flat == i)
            values = img_flat[mask_bin]

            if len(values) > 0:
                mean = np.nanmean(values)
                std = np.nanstd(values)
                err = std / np.sqrt(len(values))
            else:
                mean, err = np.nan, np.nan

            img_prof.append(mean)
            img_err.append(err)
            mod_prof.append(np.nanmean(mod_flat[mask_bin]))

        return bins, np.array(img_prof), np.array(mod_prof), np.array(img_err)

    # -------- Profile Plot --------
    def plot_profile(self, ax, r, img_prof, mod_prof, img_err):
        ax.errorbar(r, img_prof, yerr=img_err,
                    fmt='o', markersize=6, label='Image', alpha=0.7)

        ax.plot(r, mod_prof, label='Model', lw=2)

        ax.set_yscale('log')
        ax.set_xlabel("Radius (pixels)", size=15)
        ax.set_ylabel("Surface Brightness", size=15)
        ax.legend(fontsize=15)

        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)


    # -------- Histogram --------
    def plot_histogram(self, ax, sigma=3):
        res = self.residual_m[self.valid]
        res = res[~np.isnan(res)]

        mean = np.mean(res)
        std = np.std(res)
        clipped = res[np.abs(res - mean) < sigma * std]

        ax.hist(clipped, bins=50)
        ax.set_xlabel("Residual", size=15)
        ax.set_ylabel("")

        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)

    # -------- Main Plot --------
    def plot_summary(self, save_name="output.png"):
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        # Top row
        self.plot_image(axes[0, 0], self.image, "Image", cmap='viridis')
        self.plot_image(axes[0, 1], self.model, "Model", cmap='viridis')
        self.plot_residual(axes[0, 2], self.residual, "Residual")

        # Bottom row
        r, img_prof, mod_prof, img_err = self.radial_profile()
        self.plot_profile(axes[1, 0], r, img_prof, mod_prof, img_err)

        self.plot_histogram(axes[1, 1])

        axes[1, 2].axis('off')

        plt.tight_layout()
        fig.savefig(save_name)
        plt.close(fig)

if __name__=='__main__':

    plotter = FitsPlotter("O_cl1358_9.fits", "EM_cl1358_9.fits")
    plotter.plot_summary("result.png")
