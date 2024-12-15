import re
import os  # Add to existing imports at top
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                            QHBoxLayout, QLabel, QLineEdit, QComboBox, 
                            QPushButton, QTextEdit, QFrame, QFileDialog,
                            QListWidget, QAbstractItemView, QListWidgetItem)
import requests
import matplotlib.patheffects
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from astroquery.mocserver import MOCServerClass
from astroquery.hips2fits import hips2fits
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
import math

class DownloadThread(QThread):
    """
    A thread for handling asynchronous downloads of astronomical images.
    
    This class manages the download of images from various astronomical surveys using
    the HiPS (Hierarchical Progressive Survey) system. It handles coordinate resolution,
    survey querying, and image downloading without blocking the main GUI thread.
    
    Signals:
        progress (str): Emitted to report progress updates
        finished (list): Emitted when downloads complete, contains image data
        error (str): Emitted when an error occurs
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)

    def __init__(self, object_name, surveys):
        """
        Initialize the download thread.

        Parameters:
        object_name (str): The name of the astronomical object to query.
        surveys (dict): A dictionary of surveys with boolean values indicating whether to include each survey.
        """
        super().__init__()
        self.object_name = object_name
        self.surveys = surveys

    def run(self):
        """
        Execute the download thread's main operation.
        
        This method handles the complete workflow of:
        1. Resolving coordinates (from name or direct values)
        2. Querying available datasets from selected surveys
        3. Downloading images from matching datasets
        4. Emitting progress updates and final results
        """
        try:
            # Use stored coordinates instead of querying if available
            if hasattr(self, 'ra') and hasattr(self, 'dec'):
                coords = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg))
            else:
                # Original Simbad query code
                self.progress.emit(f"Querying coordinates for {self.object_name}...")
                result_table = Simbad.query_object(self.object_name)
                if result_table is None:
                    self.error.emit("Could not find object in Simbad database")
                    return
                coords = SkyCoord(result_table['RA'][0], result_table['DEC'][0], 
                                unit=(u.hourangle, u.deg))

            # Get datasets based on selected surveys
            moc_server = MOCServerClass()
            available_datasets = []

            for survey_id in self.surveys.keys():
                try:
                    self.progress.emit(f"Querying {survey_id}...")
                    # Map common survey IDs to their correct HiPS identifiers
                    survey_map = {
                        'CDS/P/2MASS/H': '2MASS/H',
                        'CDS/P/2MASS/J': '2MASS/J',
                        'CDS/P/2MASS/K': '2MASS/K',
                        'CDS/P/DSS2/red': 'DSS2/red',
                        'CDS/P/DSS2/blue': 'DSS2/blue',
                        'CDS/P/PanSTARRS/DR1/g': 'PanSTARRS/DR1/g',
                        'CDS/P/PanSTARRS/DR1/r': 'PanSTARRS/DR1/r',
                        'CDS/P/PanSTARRS/DR1/i': 'PanSTARRS/DR1/i',
                        'CDS/P/PanSTARRS/DR1/z': 'PanSTARRS/DR1/z',
                        'CDS/P/PanSTARRS/DR1/y': 'PanSTARRS/DR1/y'
                    }
                    
                    clean_id = survey_map.get(survey_id, survey_id.replace('CDS/P/', ''))
                    query = f'ID=*{clean_id}*'
                    datasets = moc_server.find_datasets(query)
                    
                    if datasets and len(datasets) > 0:
                        # Filter out multidimensional columns
                        names = [name for name in datasets.colnames if len(datasets[name].shape) <= 1]
                        df = datasets[names].to_pandas()
                        available_datasets.extend(df.to_dict('records'))
                    else:
                        self.progress.emit(f"No data found for {survey_id}, trying alternative query...")
                        # Try alternative query without wildcards
                        query = f'ID={clean_id}'
                        datasets = moc_server.find_datasets(query)
                        if datasets and len(datasets) > 0:
                            names = [name for name in datasets.colnames if len(datasets[name].shape) <= 1]
                            df = datasets[names].to_pandas()
                            available_datasets.extend(df.to_dict('records'))
                        else:
                            self.progress.emit(f"No data found for {survey_id}")
                        
                except Exception as e:
                    self.progress.emit(f"Error querying {survey_id}: {str(e)}")
                    continue

            # Download images
            images = []
            labels = []
            fov_quantity = getattr(self, 'fov', 0.1) * u.deg  # Use configured FOV or default

            for dataset in self.sort_datasets_by_filter_order(available_datasets):
                hips_id = dataset.get('ID', None)
                if hips_id is None:
                    continue

                try:
                    self.progress.emit(f"Downloading {hips_id}...")
                    img = hips2fits.query(
                        hips=hips_id,
                        ra=coords.ra,
                        dec=coords.dec,
                        fov=fov_quantity,
                        width=512,
                        height=512,
                        projection='TAN',
                        format='jpg'
                    )
                    if img is not None:
                        images.append((img, hips_id))
                except Exception as e:
                    self.progress.emit(f"Error downloading {hips_id}: {str(e)}")

            self.finished.emit(images)

        except Exception as e:
            self.error.emit(str(e))

    def sort_datasets_by_filter_order(self, datasets):
        """
        Sort datasets by their astronomical filter wavelength order.
        
        Args:
            datasets (list): List of dataset dictionaries to sort
            
        Returns:
            list: Sorted datasets from shortest to longest wavelength
        """
        # Remove duplicate filter_order dictionary and use get_filter_order_value
        return sorted(datasets, key=self.get_filter_order_value)
    
    def get_filter_order_value(self, dataset):
        """
        Get the sorting value for a dataset based on its filter.
        
        Args:
            dataset (dict): Dataset dictionary containing filter information
            
        Returns:
            int: Sorting value (lower = shorter wavelength)
        """
        filter_order = {
            'GALEX': 1,
            'PanSTARRS/DR1/g': 2,
            'PanSTARRS/DR1/r': 3,
            'PanSTARRS/DR1/i': 4,
            'PanSTARRS/DR1/z': 5,
            'PanSTARRS/DR1/y': 6,
            '2MASS/J': 7,
            '2MASS/H': 8,
            '2MASS/K': 9,
            'WISE/W1': 10,
            'WISE/W2': 11,
            'WISE/W3': 12,
            'WISE/W4': 13,
            'JWST': 14
        }
        return filter_order.get(dataset.get('ID', ''), 999)

    def download_fits(self, coords, fov_quantity, output_dir):
        """
        Download FITS format files for selected surveys.
        
        Args:
            coords (SkyCoord): Sky coordinates for the cutout center
            fov_quantity (Quantity): Field of view size in degrees
            output_dir (str): Directory path where FITS files will be saved
            
        The method creates the output directory if it doesn't exist and
        saves FITS files named according to their survey ID.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # Use the same survey mapping as in run method
        survey_map = {
            'CDS/P/2MASS/H': '2MASS/H',
            'CDS/P/2MASS/J': '2MASS/J',
            'CDS/P/2MASS/K': '2MASS/K',
            'CDS/P/DSS2/red': 'DSS2/red',
            'CDS/P/DSS2/blue': 'DSS2/blue',
            'CDS/P/PanSTARRS/DR1/g': 'PanSTARRS/DR1/g',
            'CDS/P/PanSTARRS/DR1/r': 'PanSTARRS/DR1/r',
            'CDS/P/PanSTARRS/DR1/i': 'PanSTARRS/DR1/i',
            'CDS/P/PanSTARRS/DR1/z': 'PanSTARRS/DR1/z',
            'CDS/P/PanSTARRS/DR1/y': 'PanSTARRS/DR1/y'
        }
            
        for survey_id in self.surveys.keys():
            try:
                self.progress.emit(f"Downloading FITS for {survey_id}...")
                # Clean up survey ID
                clean_id = survey_map.get(survey_id, survey_id.replace('CDS/P/', ''))
                safe_name = clean_id.replace('/', '_')
                output_file = os.path.join(output_dir, f"{safe_name}.fits")
                
                img = hips2fits.query(
                    hips=clean_id,
                    ra=coords.ra,
                    dec=coords.dec,
                    fov=fov_quantity,
                    width=512,
                    height=512,
                    projection='TAN',
                    format='fits',
                    get_query_payload=False
                )
                
                if img is not None:
                    img.writeto(output_file, overwrite=True)
                    self.progress.emit(f"Saved {output_file}")
                else:
                    self.progress.emit(f"Failed to download FITS for {survey_id}")
                    
            except Exception as e:
                self.progress.emit(f"Error downloading FITS for {survey_id}: {str(e)}")

class HipsCutoutGUI(QMainWindow):
    """
    Main GUI application for fetching and displaying astronomical image cutouts.
    
    This class provides a graphical interface for:
    - Entering object names or coordinates
    - Selecting surveys to query
    - Displaying downloaded images
    - Saving image collages
    - Downloading FITS format data
    """
    DEFAULT_OBJECT_NAME = "M 51"
    DEFAULT_CUTOUT_SIZE = "0.1"
    MOCSERVER_URL = ("http://alasky.cds.unistra.fr/MocServer/query"
                     "?expr=(hips_frame%3Dequatorial%2Cgalactic%2Cecliptic+||+"
                     "hips_frame%3D!*)+%26%26+dataproduct_type!%3Dcatalog%2Ccube+"
                     "%26%26+hips_service_url%3D*&get=record")
    DEFAULT_SURVEYS = [
        "CDS/P/2MASS/color",
        "CDS/P/HST/EPO",
        "CDS/P/SDSS9/color"
    ]

    def __init__(self):
        super().__init__()
        self.setWindowTitle("HIPS Cutout Viewer")
        self.setGeometry(100, 100, 1200, 800)
        self.selected_surveys = []
        self.hips_data = self.fetch_hips_data()
        self.init_ui()

    def fetch_hips_data(self):
        """
        Fetch available HiPS surveys from the MOCServer.
        
        Returns:
            dict: Dictionary of survey IDs and their display names
            
        The method queries the MOCServer for available surveys and parses
        the response to create a mapping of survey identifiers.
        """
        try:
            response = requests.get(self.MOCSERVER_URL)
            if response.status_code != 200:
                print(f"Error fetching data: {response.status_code}")
                return {}

            # Parse the response text for lines containing "ID = "
            survey_ids = {}
            for line in response.text.split('\n'):
                if 'ID                   = ' in line:
                    survey_id = line.split('= ')[1].strip()
                    if survey_id:  # Only add non-empty IDs
                        survey_ids[survey_id] = survey_id  # Use full ID as display name

            return survey_ids

        except Exception as e:
            print(f"Error fetching HIPS data: {e}")
            return {}

    def init_ui(self):
        self.create_central_widget()
        self.create_input_area()
        self.create_progress_log()
        self.create_plot_area()
        self.download_thread = None
        
        # Add default surveys to selected list
        default_surveys = [
            "CDS/P/2MASS/color",
            "CDS/P/HST/EPO",
            "CDS/P/SDSS9/color"
        ]
        
        for survey_id in default_surveys:
            if survey_id in self.hips_data:
                item = QListWidgetItem(survey_id)
                item.setData(Qt.ItemDataRole.UserRole, survey_id)
                self.selected_list.addItem(item)
                self.selected_surveys.append(survey_id)

    def create_central_widget(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.layout = QVBoxLayout(central_widget)

    def create_input_area(self):
        input_widget = QWidget()
        input_layout = QVBoxLayout(input_widget)
        inputs_layout = self.create_inputs_layout()
        input_layout.addLayout(inputs_layout)
        survey_layout = self.create_survey_selection()  # Replace checkbox_layout
        input_layout.addLayout(survey_layout)
        button_layout = self.create_button_layout()
        input_layout.addLayout(button_layout)
        self.layout.addWidget(input_widget)

    def create_inputs_layout(self):
        inputs_layout = QHBoxLayout()
        
        # Object name section
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Object Name:"))
        self.object_name = QLineEdit(self.DEFAULT_OBJECT_NAME)
        self.object_name.setFixedWidth(200)
        self.object_name.textChanged.connect(self.clear_coords)
        name_layout.addWidget(self.object_name)
        inputs_layout.addLayout(name_layout)
        inputs_layout.addSpacing(20)

        # RA section
        self.ra_input = QLineEdit()
        self.ra_input.setFixedWidth(100)
        self.ra_input.setPlaceholderText("RA (deg)")
        inputs_layout.addWidget(QLabel("RA:"))
        inputs_layout.addWidget(self.ra_input)

        # Dec section
        self.dec_input = QLineEdit()
        self.dec_input.setFixedWidth(100)
        self.dec_input.setPlaceholderText("Dec (deg)")
        inputs_layout.addWidget(QLabel("Dec:"))
        inputs_layout.addWidget(self.dec_input)

        # Add cutout size section
        inputs_layout.addSpacing(20)
        inputs_layout.addWidget(QLabel("Size (deg):"))
        self.size_input = QLineEdit(self.DEFAULT_CUTOUT_SIZE)
        self.size_input.setFixedWidth(60)
        inputs_layout.addWidget(self.size_input)

        # Resolve button
        self.resolve_btn = QPushButton("Resolve Name")
        self.resolve_btn.clicked.connect(self.resolve_coordinates)
        inputs_layout.addWidget(self.resolve_btn)
        
        inputs_layout.addStretch()
        return inputs_layout

    def create_survey_selection(self):
        """
        Create the survey selection interface components.
        
        Returns:
            QHBoxLayout: Layout containing the survey selection widgets
            
        Creates a two-column layout with available surveys on the left
        and selected surveys on the right, with add/remove buttons.
        """
        survey_layout = QHBoxLayout()
        
        # Left side - Dropdown and Add button
        left_layout = QVBoxLayout()
        self.survey_combo = QComboBox()
        # Add items sorted by full ID
        sorted_items = sorted(self.hips_data.items())
        for survey_id, _ in sorted_items:
            self.survey_combo.addItem(survey_id, survey_id)  # Both display and data are full ID
        left_layout.addWidget(QLabel("Available Surveys:"))
        left_layout.addWidget(self.survey_combo)
        
        add_btn = QPushButton("Add →")
        add_btn.clicked.connect(self.add_survey)
        left_layout.addWidget(add_btn)
        left_layout.addStretch()
        
        # Right side - Selected surveys list
        right_layout = QVBoxLayout()
        right_layout.addWidget(QLabel("Selected Surveys:"))
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        right_layout.addWidget(self.selected_list)
        
        remove_btn = QPushButton("← Remove")
        remove_btn.clicked.connect(self.remove_survey)
        right_layout.addWidget(remove_btn)
        
        survey_layout.addLayout(left_layout)
        survey_layout.addLayout(right_layout)
        return survey_layout

    def add_survey(self):
        """
        Add the currently selected survey to the active surveys list.
        
        Checks for duplicates and updates the selected_surveys list
        with the newly added survey's ID.
        """
        current_index = self.survey_combo.currentIndex()
        display_name = self.survey_combo.currentText()
        survey_id = self.survey_combo.itemData(current_index)
        
        if display_name not in [self.selected_list.item(i).text() 
                              for i in range(self.selected_list.count())]:
            item = QListWidgetItem(display_name)
            item.setData(Qt.ItemDataRole.UserRole, survey_id)
            self.selected_list.addItem(item)
            self.selected_surveys = [self.selected_list.item(i).text() 
                                   for i in range(self.selected_list.count())]

    def remove_survey(self):
        """
        Remove the selected survey from the active surveys list.
        
        Updates the selected_surveys list to reflect the removal.
        """
        current_row = self.selected_list.currentRow()
        if current_row >= 0:
            self.selected_list.takeItem(current_row)
            self.selected_surveys = [self.selected_list.item(i).text() 
                                   for i in range(self.selected_list.count())]

    def create_button_layout(self):
        button_layout = QHBoxLayout()
        self.submit_btn = QPushButton("Get Cutouts")
        self.submit_btn.clicked.connect(self.submit)
        button_layout.addWidget(self.submit_btn)
        self.save_btn = QPushButton("Save Collage")
        self.save_btn.clicked.connect(self.save_collage)
        self.save_btn.setEnabled(False)
        button_layout.addWidget(self.save_btn)
        self.reset_btn = QPushButton("Reset")
        self.reset_btn.clicked.connect(self.reset)
        button_layout.addWidget(self.reset_btn)
        
        self.fits_btn = QPushButton("Download FITS")
        self.fits_btn.clicked.connect(self.download_fits)
        self.fits_btn.setEnabled(False)
        button_layout.addWidget(self.fits_btn)
        
        return button_layout

    def create_progress_log(self):
        self.log_text = QTextEdit()
        self.log_text.setMaximumHeight(100)
        self.log_text.setReadOnly(True)
        self.layout.addWidget(self.log_text)

    def create_plot_area(self):
        """Create the plot area for displaying images"""
        self.plot_widget = QFrame()
        self.layout.addWidget(self.plot_widget)
        self.plot_layout = QVBoxLayout(self.plot_widget)

    def log(self, message):
        """
        Append a message to the log text area.
        
        Parameters:
        message (str): The message to append to the log.
        """
        self.log_text.append(message)

    def clear_coords(self):
        """Clear coordinates when object name is modified"""
        self.ra_input.clear()
        self.dec_input.clear()

    def resolve_coordinates(self):
        """
        Resolve the object name to its coordinates using the Simbad database.
        """
    def submit(self):
        """Modified submit to handle both name and coordinate-based queries"""
        self.submit_btn.setEnabled(False)
        self.log_text.clear()

        # Check if we have direct coordinates
        try:
            if self.ra_input.text() and self.dec_input.text():
                ra = float(self.ra_input.text())
                dec = float(self.dec_input.text())
            else:
                # If no coordinates, try to resolve name first
                self.resolve_coordinates()
                if not self.ra_input.text() or not self.dec_input.text():  # Fixed extra parenthesis
                    self.submit_btn.setEnabled(True)
                    return
                ra = float(self.ra_input.text())
                dec = float(self.dec_input.text())

            # Get selected surveys using stored survey IDs
            selected_surveys = {}
            for i in range(self.selected_list.count()):
                item = self.selected_list.item(i)
                survey_id = item.data(Qt.ItemDataRole.UserRole)
                if survey_id:
                    selected_surveys[survey_id] = True

            if not selected_surveys:
                self.log("Please select at least one survey")
                self.submit_btn.setEnabled(True)
                return

            # Get cutout size
            try:
                fov = float(self.size_input.text())
            except ValueError:
                self.log("Invalid cutout size. Using default value of 0.1 degrees.")
                fov = 0.1

            # Start download thread with coordinates and FOV
            self.download_thread = DownloadThread(f"RA: {ra}, Dec: {dec}", selected_surveys)
            self.download_thread.ra = ra
            self.download_thread.dec = dec
            self.download_thread.fov = fov
            self.download_thread.progress.connect(self.log)
            self.download_thread.finished.connect(self.display_images)
            self.download_thread.error.connect(self.handle_error)
            self.download_thread.start()

        except ValueError:
            self.log("Invalid coordinate format. Please enter valid decimal degrees.")
            self.submit_btn.setEnabled(True)

    def handle_error(self, error_message):
        """
        Handle errors by logging the error message and disabling the submit and save buttons.

        Parameters:
        error_message (str): The error message to log.
        """
        self.log(f"Error: {error_message}")
        self.submit_btn.setEnabled(True)
        self.save_btn.setEnabled(False)  # Disable save button on error

    def resolve_coordinates(self):
        """Resolve object name to coordinates"""
        
        try:
            result_table = Simbad.query_object(self.object_name.text())
            if result_table is None:
                self.log("Could not find object in Simbad database")
                return
            
            coords = SkyCoord(result_table['RA'][0], result_table['DEC'][0], 
                            unit=(u.hourangle, u.deg))
            
            # Update coordinate fields
            self.ra_input.setText(f"{coords.ra.deg:.6f}")
            self.dec_input.setText(f"{coords.dec.deg:.6f}")
            self.log(f"Resolved coordinates: RA={coords.ra.deg:.6f}°, Dec={coords.dec.deg:.6f}°")
            
        except Exception as e:
            self.log(f"Error resolving coordinates: {str(e)}")

    def submit(self):
        """Modified submit to handle both name and coordinate-based queries"""
        self.submit_btn.setEnabled(False)
        self.log_text.clear()

        # Check if we have direct coordinates
        try:
            if self.ra_input.text() and self.dec_input.text():
                ra = float(self.ra_input.text())
                dec = float(self.dec_input.text())
            else:
                # If no coordinates, try to resolve name first
                self.resolve_coordinates()
                if not self.ra_input.text() or not self.dec_input.text():
                    self.submit_btn.setEnabled(True)
                    return
                ra = float(self.ra_input.text())
                dec = float(self.dec_input.text())

            # Get selected surveys
            selected_surveys = {survey: True for survey in self.selected_surveys}
            
            if not selected_surveys:
                self.log("Please select at least one survey")
                self.submit_btn.setEnabled(True)
                return

            # Get cutout size
            try:
                fov = float(self.size_input.text())
            except ValueError:
                self.log("Invalid cutout size. Using default value of 0.1 degrees.")
                fov = 0.1

            # Start download thread with coordinates and FOV
            self.download_thread = DownloadThread(f"RA: {ra}, Dec: {dec}", selected_surveys)
            self.download_thread.ra = ra
            self.download_thread.dec = dec
            self.download_thread.fov = fov
            self.download_thread.progress.connect(self.log)
            self.download_thread.finished.connect(self.display_images)
            self.download_thread.error.connect(self.handle_error)
            self.download_thread.start()

        except ValueError:
            self.log("Invalid coordinate format. Please enter valid decimal degrees.")
            self.submit_btn.setEnabled(True)

    def handle_error(self, error_message):
        """
        Handle errors by logging the error message and disabling the submit and save buttons.

        Parameters:
        error_message (str): The error message to log.
        """
        self.log(f"Error: {error_message}")
        self.submit_btn.setEnabled(True)
        self.save_btn.setEnabled(False)  # Disable save button on error

    def display_images(self, image_data):
        # Clear existing plot
        for i in reversed(range(self.plot_layout.count())): 
            self.plot_layout.itemAt(i).widget().setParent(None)

        # Filter out any None or empty images
        valid_images = [(img, label) for img, label in image_data if img is not None]
        
        if not valid_images:
            self.log("No images were downloaded.")
            self.submit_btn.setEnabled(True)
            return

        # Create plot
        num_images = len(valid_images)
        cols = 5
        rows = math.ceil(num_images / cols)

        fig = Figure(figsize=(15, 3 * rows))
        for i, (img, label) in enumerate(valid_images):
            ax = fig.add_subplot(rows, cols, i+1)
            ax.imshow(img, origin='lower', cmap='gray')
            ax.set_title(label, fontsize=8)
            self.add_markers(ax, img)
            ax.axis('off')

        fig.tight_layout()

        canvas = FigureCanvasQTAgg(fig)
        self.plot_layout.addWidget(canvas)
        self.submit_btn.setEnabled(True)
        
        # Store figure reference and enable save button
        self.current_figure = fig
        self.save_btn.setEnabled(True)
        self.fits_btn.setEnabled(True)  # Enable FITS button after successful image display

    def save_collage(self):
        if hasattr(self, 'current_figure'):
            filename, _ = QFileDialog.getSaveFileName(
                self,
                "Save Collage",
                f"{self.object_name.text().replace(' ', '_')}_collage.jpg",
                "JPEG Files (*.jpg *.jpeg)"
            )
            if filename:
                self.log(f"Saving collage to {filename}")
                try:
                    self.current_figure.savefig(filename, 
                                              format='jpg', 
                                              bbox_inches='tight', 
                                              dpi=300)
                    self.log("Collage saved successfully!")
                except Exception as e:
                    self.log(f"Error saving collage: {str(e)}")

    def add_markers(self, ax, img):
        """
        Add orientation markers and scale bar to an image plot.
        
        Args:
            ax (Axes): Matplotlib axes object to add markers to
            img (array): Image data used for dimensioning markers
            
        Adds:
        - North-East orientation arrows
        - Scale bar showing 1 arcminute
        - Labels for orientation and scale
        """
        # Add North-East arrows
        arrow_length = int(0.1 * img.shape[1])
        x_base = int(0.15 * img.shape[1])
        y_base = int(0.15 * img.shape[0])
        
        # North-East arrows
        ax.arrow(x_base, y_base, 0, arrow_length, color='white', 
                head_width=5, head_length=10, fc='white', ec='white')
        ax.text(x_base, y_base + arrow_length + 15, 'N', color='white', 
                horizontalalignment='center')
        
        ax.arrow(x_base, y_base, arrow_length, 0, color='white', 
                head_width=5, head_length=10, fc='white', ec='white')
        ax.text(x_base + arrow_length + 15, y_base, 'E', color='white', 
                verticalalignment='center')
        
        # Add scale bar (1 arcmin)
        scale_pixels = int(512 / 6)  # 1 arcmin (1/6 of degree) in pixels
        bar_y = int(0.1 * img.shape[0])  # 10% from bottom
        bar_x = int(0.9 * img.shape[1])  # 90% from left
        
        # Draw white scale bar with black edge for better visibility
        ax.plot([bar_x - scale_pixels, bar_x], [bar_y, bar_y], 
                'k-', linewidth=4)  # black edge
        ax.plot([bar_x - scale_pixels, bar_x], [bar_y, bar_y], 
                'w-', linewidth=2)  # white center
        
        # Add scale text with black edge for better visibility
        ax.text(bar_x - scale_pixels/2, bar_y + 15, "1'", 
                color='white', fontsize=8,
                horizontalalignment='center',
                verticalalignment='bottom',
                path_effects=[
                    matplotlib.patheffects.withStroke(linewidth=2, 
                                                    foreground='black')
                ])

    def reset(self):
        """Reset inputs to default values while preserving default surveys"""
        self.object_name.setText(self.DEFAULT_OBJECT_NAME)
        self.ra_input.clear()
        self.dec_input.clear()
        self.size_input.setText(self.DEFAULT_CUTOUT_SIZE)
        
        # Clear and restore default surveys
        self.selected_list.clear()
        self.selected_surveys = []
        
        # Re-add default surveys
        for survey_id in self.DEFAULT_SURVEYS:
            if survey_id in self.hips_data:
                item = QListWidgetItem(survey_id)
                item.setData(Qt.ItemDataRole.UserRole, survey_id)
                self.selected_list.addItem(item)
                self.selected_surveys.append(survey_id)
        
        # Clear the log and any displayed images
        self.log_text.clear()
        for i in reversed(range(self.plot_layout.count())): 
            self.plot_layout.itemAt(i).widget().setParent(None)
        
        self.save_btn.setEnabled(False)
        self.submit_btn.setEnabled(True)
        self.fits_btn.setEnabled(False)
        
        self.log("All inputs reset to defaults")

    def download_fits(self):
        """
        Initiate FITS format file downloads for selected surveys.
        
        Validates inputs, prepares download parameters, and starts
        a download thread to fetch FITS files for each selected survey.
        Files are saved in a 'fits' subdirectory of the current working directory.
        """
        try:
            if not (self.ra_input.text() and self.dec_input.text()):
                self.log("Please resolve coordinates first")
                return
                
            ra = float(self.ra_input.text())
            dec = float(self.dec_input.text())
            
            # Get selected surveys
            selected_surveys = {}
            for i in range(self.selected_list.count()):
                item = self.selected_list.item(i)
                survey_id = item.data(Qt.ItemDataRole.UserRole)
                if survey_id:
                    selected_surveys[survey_id] = True
            
            if not selected_surveys:
                self.log("Please select at least one survey")
                return

            # Get cutout size
            try:
                fov = float(self.size_input.text())
            except ValueError:
                self.log("Invalid cutout size. Using default value of 0.1 degrees.")
                fov = 0.1

            # Start download thread
            self.fits_btn.setEnabled(False)
            self.download_thread = DownloadThread(f"RA: {ra}, Dec: {dec}", selected_surveys)
            self.download_thread.ra = ra
            self.download_thread.dec = dec
            self.download_thread.fov = fov
            
            # Create fits directory in current working directory
            fits_dir = os.path.join(os.getcwd(), 'fits')
            
            self.download_thread.progress.connect(self.log)
            self.download_thread.finished.connect(lambda: self.fits_download_complete())
            self.download_thread.error.connect(self.handle_error)
            
            # Call download_fits method
            self.download_thread.download_fits(
                SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg)),
                fov * u.deg,
                fits_dir
            )
            
        except Exception as e:
            self.log(f"Error starting FITS download: {str(e)}")
            self.fits_btn.setEnabled(True)

    def fits_download_complete(self):
        """
        Handle the completion of FITS downloads.
        
        Re-enables the FITS download button and logs completion message.
        """
        self.fits_btn.setEnabled(True)
        self.log("FITS downloads complete")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = HipsCutoutGUI()
    window.show()
    sys.exit(app.exec())
