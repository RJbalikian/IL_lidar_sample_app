import pathlib
import pyproj
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray as rxr
import streamlit as st
from rasterio import MemoryFile
from owslib.wms import WebMapService
from io import BytesIO
import plotly.graph_objects as go
import plotly.express as px

#matplotlib.use("agg")

CRS_LIST = pyproj.database.query_crs_info()
CRS_STR_LIST = [f"{crs.auth_name}:{crs.code} - {crs.name}" for crs in CRS_LIST]
CRS_DICT = {f"{crs.auth_name}:{crs.code} - {crs.name}": crs for crs in CRS_LIST}
LIDAR_SERVICE_URL=r"https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WMSServer?request=GetCapabilities&service=WMS"
RASTER_SRC_DICT = {"ISGS Statewide Lidar":LIDAR_SERVICE_URL,
                   "Other Web Service":"get_service_info",
                   "Raster file":"get_file_name"
                   }

DEFAULT_POINTS_CRS = "EPSG:6345 - NAD83(2011) / UTM zone 16N"
DEFAULT_POINTS_CRS_INDEX = CRS_STR_LIST.index(DEFAULT_POINTS_CRS)

DEFAULT_OUTPUT_CRS = DEFAULT_POINTS_CRS

def main():
    print("Start Main Loop")
    #st.container(key='main_container')

    with st.container():
        if hasattr(st.session_state, 'elev_fig'):
            st.plotly_chart(st.session_state.elev_fig)
            st.session_state.coords_df
            
    with st.sidebar:
        print('Building sidebar')
        st.title("Raster Sampling on the Web")

        st.selectbox(label="Select raster/elevation data source",
                     options=["ISGS Statewide Lidar", "Other Web Service", "Raster file"],
                     index=0, key='raster_source_select',
                     help="Select the source you would like to use for elevation.",
                     on_change=on_raster_source_change,
                     disabled=True,
                     )

        st.selectbox(label="Raster CRS",
                     options=CRS_STR_LIST,
                     index=CRS_STR_LIST.index("EPSG:3857 - WGS 84 / Pseudo-Mercator"),
                     key='raster_crs',
                     disabled=True)


        #raster_expander = st.expander(label='Raster Info.')
        
        st.segmented_control(label='Select point type',
                             options=["Enter coords.",
                                      "CSV File"],
                             selection_mode='single',
                             default="Enter coords.",
                             key="points_source",
                             )

        xcoordCol, ycoordCol = st.columns(spec=[0.5, 0.5], gap='small', border=False)
    
        xcoordCol.text_input(label='X Coordinate',
                      on_change=on_xcoord_change,
                      placeholder="X Coord Value",
                      help="If you enter two comma or space/tab separated values, it will automatically parse the 2nd as your ycoord.",
                      key='xcoord'
                      )

        ycoordCol.text_input(label='Y Coordinate',
                      placeholder="Y Coord Value",
                      key='ycoord'
                      )        
        
        st.selectbox(label="Point CRS",
                     options=CRS_STR_LIST,
                     index=DEFAULT_POINTS_CRS_INDEX,
                     key='point_crs')

        #point_expander = st.expander(label='Point Info.')

        st.divider()
        st.selectbox(label="Output CRS",
                     options=CRS_STR_LIST,
                     index=DEFAULT_POINTS_CRS_INDEX,
                     key='output_crs')

        st.button(label="Sample Elevation Data",
                  key="sample_data_button",
                  icon=":material/background_dot_small:",
                  on_click=get_elevation,
                  type='primary',
                  )
        print('Done Building sidebar')
    print("MAIN LOOP END")


# FUNCTIONS
def on_raster_source_change():
    return

def on_xcoord_change():
    sepList = [',', ' ', '\t']
    if hasattr(st.session_state, 'xcoord'):
        st.session_state.xcoord = st.session_state.xcoord.strip()
        for sep in sepList:
            if sep in st.session_state.xcoord:
                xCoordList = str(st.session_state.xcoord).split(sep)

                if len(xCoordList) > 1:
                    st.session_state.xcoord = xCoordList[0].strip()
                    st.session_state.ycoord = xCoordList[1].strip()
                break

def get_elevation(coords=None, 
                  elevation_col_name='elevation', 
                  xcoord_col_name='xcoord', ycoord_col_name='ycoord', 
                  points_crs=None, output_crs=None,
                  elevation_source=None, elev_source_type='service', 
                  raster_crs=None, show_plot=True):

    print("GETTING ELEVATION")

    if coords is None:
        if st.session_state.points_source=='Enter coords.':
            coords = (st.session_state.xcoord, st.session_state.ycoord)
            #coords = (-88.857362, 42.25637743)

    # Leave everything below here alone
    if elevation_source is None:
        elevation_source = LIDAR_SERVICE_URL#RASTER_SRC_DICT[st.session_state['raster_source_select']]

    if points_crs is None:
        points_crs = int(CRS_DICT[st.session_state.point_crs].code)
        points_crs_name = CRS_DICT[st.session_state.point_crs].name

    if raster_crs is None:
        raster_crs = int(CRS_DICT[st.session_state.raster_crs].code)
        raster_crs_name = CRS_DICT[st.session_state.raster_crs].name

    if output_crs is None:
        output_crs = int(CRS_DICT[st.session_state.output_crs].code)
        output_crs_name = CRS_DICT[st.session_state.output_crs].name
    
    elev_wms = wms_statewide_lidar_url = elevation_source

    ptCoordTransformerOUT = pyproj.Transformer.from_crs(crs_from=points_crs, crs_to=output_crs, always_xy=True)
    ptCoordTransformerRaster = pyproj.Transformer.from_crs(crs_from=points_crs, crs_to=raster_crs, always_xy=True)
    wms_statewide_lidar_url = elev_wms

    if isinstance(coords, (tuple, list)):
        xcoord, ycoord = coords

        xcoord_OUT, ycoord_OUT = ptCoordTransformerOUT.transform(xcoord, ycoord)
        xcoord_RAST, ycoord_RAST = ptCoordTransformerRaster.transform(xcoord, ycoord)

        minX = maxX = xcoord_OUT
        minY = maxY = ycoord_OUT

        minXRast = maxXRast= xcoord_RAST
        minYRast = maxYRast = ycoord_RAST

        cols = [f"{points_crs}_xIN", f"{points_crs}_yIN", f"{output_crs}_x", f"{output_crs}_y"]
        coords = pd.DataFrame([[xcoord, ycoord, xcoord_OUT, ycoord_OUT]], columns=cols)

    elif isinstance(coords, (pd.DataFrame, gpd.GeoDataFrame, str, pathlib.Path)):
        if isinstance(coords, (pathlib.Path, str)):
            coords = pd.read_csv(coords)
        xcoord = coords[xcoord_col_name]
        ycoord = coords[ycoord_col_name]

        if str(raster_crs) == '3857':
            tempCoord = ycoord
            ycoord = xcoord
            xcoord = tempCoord

        xcoord_OUT, ycoord_OUT = ptCoordTransformerOUT.transform(xcoord, ycoord)

        minX = min(xcoord_OUT)
        minY = min(ycoord_OUT)
        maxX = max(xcoord_OUT)
        maxY = max(ycoord_OUT)

    xPad = (maxXRast-minXRast)*0.1
    yPad = (maxYRast-minYRast)*0.1

    if float(xPad) == 0.0:
        xPad = maxXRast * 0.01
        xPad = abs(xPad)
        if abs(xPad) > 1000:
            xPad = 1000

    if float(yPad) == 0.0:
        yPad = maxYRast * 0.01
        yPad = abs(yPad)
        if yPad > 1000:
            yPad = 1000

    if elev_source_type=='service':
        wms = WebMapService(wms_statewide_lidar_url)

        layer_name = '0'
        layer = wms[layer_name]

        bbox = (minXRast-xPad, minYRast-yPad, maxXRast+xPad, maxYRast+yPad)
        img = wms.getmap(
            layers=['IL_Statewide_Lidar_DEM_WGS:None'],
            srs='EPSG:3857',
            bbox=bbox,
            size=(256, 256),
            format='image/tiff',
            transparent=True
            )

        bio = BytesIO(img.read())
        lidarData_rxr = rxr.open_rasterio(bio)
        lidarData = lidarData_rxr.rio.reproject(output_crs)
        if 'band' in lidarData.dims:
            lidarData = lidarData.isel(band=0)

    elif elev_source_type == 'file':
        lidarData_rxr = rxr.open_rasterio(elevation_source)
        lidarData = lidarData_rxr.sel(x=slice(minX-xPad, maxX+xPad), y=slice(maxY+yPad, minY-yPad))
        lidarData = lidarData.rio.reproject(output_crs)
    
    minLidarVal = lidarData.min().values
    maxLidarVal = lidarData.max().values
    lidarValRange = maxLidarVal - minLidarVal
    
    vMin = minLidarVal + 0.2*lidarValRange
    vMax = maxLidarVal - 0.2*lidarValRange

    # Elevation in meters
    lidarData_m = lidarData * 0.3048

    print(lidarData_m)
    # Calculate elevation and add to df
    elev_m = []
    for i, row in coords.iterrows():
        lidarSel = lidarData_m.copy()
        xc = row[f"{output_crs}_x"]
        yc = row[f"{output_crs}_y"]

        elevVal = lidarSel.sel(x=xc, y=yc, 
                               method='nearest').values

        elev_m.append(elevVal)
        
    coords['Elev_m'] = elev_m
    coords['Elev_ft'] = coords['Elev_m'] / 0.3048

    if show_plot:
        lidarDA = lidarData
        data = lidarDA.values
        x_coords = lidarDA.x.values  
        y_coords = lidarDA.y.values

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=data,
            x=x_coords,
            y=y_coords,
            colorscale='Geyser',
            zmin=vMin,
            zmax=vMax,
            name='Elevation',
        ))
                
        #[print(f"{c['Elev_ft']} ft<br>{c['Elev_m']} m") for i, c in coords.iterrows()]
        # Add the point marker
        fig.add_trace(go.Scatter(
            x=coords[f"{output_crs}_x"],
            y=coords[f"{output_crs}_y"],
            mode='markers+text',
            marker=dict(
                symbol='star',
                size=15,
                color='red',
                line=dict(width=2, color='black')
            ),
            text=[f"{float(c['Elev_ft']):.1f} ft<br>{float(c['Elev_m']):.2f} m" for i, c in coords.iterrows()],
            textposition="top right",
            textfont=dict(color="red", size=12),
            name=f"Point Elevation",
            #showlegend=True
        ))

        if len(str(int(minX))) > 4:
            tickFormat = ".0f"
        else:
            tickFormat = ".2f"

        # Update layout
        fig.update_layout(
            title="Elevation Data",
            xaxis_title="X Coordinate",
            yaxis_title="Y Coordinate",
            yaxis_tickformat=tickFormat,
            xaxis_tickformat=tickFormat,
            width=800,
            height=800,
            xaxis=dict(scaleanchor='y'),
            autosize=True,
            #showlegend=True,
            coloraxis_colorbar=dict(#yanchor="top", y=1, x=0,
                                    orientation='h',
                                    ticks="outside",
        ))

        #fig.show()
        #st.session_state.main_container.pyplot(fig)
        #with st.session_state.main_container:
        #    st.plotly_chart(fig)

    st.session_state.elev_fig = fig
    st.session_state.coords_df = coords
    st.balloons()
    return #fig, coords


if __name__ == "__main__":
    main()